import os
import sys
import numpy as np 
import pandas as pd
import itertools
from itertools import islice,product,starmap
import multiprocessing as mp
from collections import OrderedDict
import re
import random
from random import sample
import time
import networkx as nx
import bz2
import json
from multiprocess import Process, Manager
from scipy.stats import fisher_exact, hmean
import scipy
import csv
from Bio import Entrez
from Bio.Entrez import efetch,esearch, read
from os import listdir
from os.path import isfile, join
from collections import Counter
import glob
from operator import itemgetter
try:
	import _pickle as pickle
except ModuleNotFoundError:
	import pickle


# Pickle a file and then compress it into a file with extension 
def compress_pickle(title, data):
	with bz2.BZ2File(title + ".pbz2", "w") as f: 
		pickle.dump(data, f)


#
# decompress .pbz2 file to a variable
def decompress_pickle(file):
	data = bz2.BZ2File(file, "rb") 
	data = pickle.load(data)
	return data


#
print("modules imported")
#
wd= sys.argv[1]
#
tag= sys.argv[2] #iter1
#
sd=re.split("tmp",wd)[0]
#
tmp=''.join(re.split("\/"+tag,wd))
#
#
prev_id_file=tmp+"source_pmids.txt"
comb=int(sys.argv[3])
P=int(sys.argv[4])
#
it=sys.argv[5]
#
# input resources
#aliasdf=pd.read_csv(sd+"input/gene_aliases_df.tsv",sep='\t')
aliasdf=pd.read_csv(sd+'input/aliases_uniprotensembl_as_df_Paralogues.tsv',sep='\t')
# remove bad aliases
aliases=aliasdf.alias_symbol.tolist()
# unless they are reference genes => regardless being a reference gene
bad_hits=list(set(filter(lambda x: re.match("^[a-zA-Z]{1,2}$",str(x)),aliases)))
aliasdf=aliasdf[~aliasdf['alias_symbol'].isin(bad_hits)]
#defining functions
Entrez.email='Luca.musella@uk-erlangen.de'
Entrez.api_key='b89451c47164b8705c76108eec6386658607'
# GET #RETMAX PMIDS BY MESH TERM
mesh_dict=decompress_pickle(sd+'input/meshD_names_to_ids.pbz2')
ref_mesh={key.replace(',','').lower():key for key in mesh_dict.keys()}
#
# load pmids that were found before
with open(prev_id_file,'r') as f:
	prev_ids=f.readlines()


prev_ids=[s for s in prev_ids if not s.startswith('#')]
prev_ids=[s for s in prev_ids if not s.startswith('\n')]
prev_ids=[int(s.split('\n')[0]) for s in prev_ids]
prev_ids=list(set(prev_ids))
#
try:
	searched=decompress_pickle(sd+'input/searched_queries.pbz2')
except:
	searched=pd.Series(dtype=object)
#
#print(searched)
# load the "complete graph"
comp_node_tbl=pd.read_csv(wd+tag+"_Complete_nodes_table_subgraph.tsv",delimiter='\t')
comp_edge_tbl=pd.read_csv(wd+tag+"_Complete_edges_table_subgraph.tsv",delimiter='\t')
#
cols=list(comp_node_tbl.columns)
# category/label is 2nd column
cat=cols[1]
#
attr={comp_node_tbl[cols[0]][i]:{cols[j]:comp_node_tbl[cols[j]][i] for j in range(1,len(cols))} for i in range(len(comp_node_tbl))}
#
comp_Net=nx.from_pandas_edgelist(comp_edge_tbl,'kw1','kw2',['pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration'])
# adding node attributes
nx.set_node_attributes(comp_Net,attr)
del attr
#
ents=comp_node_tbl['KW'].tolist() # all entities (genes and MeSH)
#
print("Complete Network and attributes loaded")
#
# dictionary of sub-gene-graphs (1 for edges 1 for nodes)
nod_fil=glob.glob(wd+'data/gene-subgraphs/*nodes*')
edg_fil=glob.glob(wd+'data/gene-subgraphs/*edges*')
#
nod_fil.sort()
edg_fil.sort()
#
#
subg_noded={i:pd.read_csv(j,delimiter='\t') for i,j in zip(range(0,len(nod_fil)),nod_fil)}
subg_edged={i:pd.read_csv(j,delimiter='\t') for i,j in zip(range(0,len(edg_fil)),edg_fil)}
#
def recnet(nod,edg):
	# nod and edg are dict iterms of the form (key,df)
	key,noddf =nod #list(subg_noded.items())[0] #nod
	key,edgdf =edg #list(subg_edged.items())[0] #edg
	cols=list(noddf.columns)
	# category/label is 2nd column
	cat=cols[1]
	#
	attr={noddf[cols[0]][i]:{cols[j]:noddf[cols[j]][i] for j in range(1,len(cols))} for i in range(len(noddf))}
	#
	#print(attr)
	#
	net=nx.from_pandas_edgelist(edgdf,'kw1','kw2',['pg','crit_p','OR','edge_weight','inv_edge_weight','iteration'])
	# adding node attributes
	nx.set_node_attributes(net,attr)
	#print(net.nodes.data())
	del attr
	#
	return (key,net)
# #
# k,G=recnet(list(subg_noded.items())[0],list(subg_edged.items())[0])
# subgnets={}
# for i in range(0,len(nod_fil)):
# 	k,G=recnet(list(subg_noded.items())[i],list(subg_edged.items())[i])
# 	subgnets[k]=G
# #
subgnets=dict(starmap(recnet,zip(subg_noded.items(),subg_edged.items())))
# REMOVE SUB-GENE GRAPHS IF ABSENT FROM THE COMPLETE NETWORK (NO PATH CAN BE FOUND)
subgnets={k:v for k,v in subgnets.items() if len(v.nodes & comp_Net.nodes) > 0}
#
### here comes the convoluted part ###
print("subgraphs have been generated")
# STRATEGY #
#
# 1) function that applies to 2 graphs at a time
# 2) create a sorted list of node pairs, ranked by respective average rank
# 3) search for shortest path, passing through MeSH (no need to check, since
#	 the nodes belong to 2 subraphs)
# 4) repeat for any combination 
#
# inter_query, adapted from "query_motif_mix.py"
#
def inter_query(path,searched,comb,init=1,ref_mesh=ref_mesh,ref_gen=aliasdf,prev_ids=prev_ids):
	#global searched
	# terms is a concatenation of two motifs
	db="pubmed"
	retmax=95000 # max allowed 100k
	idtype="acc"
	retstart=init
	query=str()
	def esearch_handle(db,retstart,retmax,term,idtype):
		res=[]
		count=2
		# catch papers above the retmax limit
		while count>retstart:
			try:
				print("search for %s: " % (term),count, retstart)
				#print("search for %s",count,retstart) % (term)
				handle = Entrez.esearch(db=db,retstart=retstart, retmax=retmax, term=term, idtype=idtype)
				record = Entrez.read(handle)
				handle.close()
			except:
				print("Internet connection problems?")
				time.sleep(5)
			res += record['IdList']
			count=int(record['Count']) # actually immutable
			retstart=1+len(res)	
		#
		res=[int(idd) for idd in res]				       
		return res
	# updated searched series if necessary
	for i in path:
		if i not in searched.index:
			search=dict({i:[]}) # then new search 
			if i in aliasdf.symbol.values: # it's a gene							
				# detect then all aliases
				aliases=aliasdf[aliasdf['symbol']==i]['alias_symbol'].tolist()
				aliases=['+'.join(al.split())+"[Title/Abstract] NOT review[Publication Type]" for al in aliases]
				j=0
				while j < len(aliases):
					# fill spaces 
					search[i]+=esearch_handle(db,retstart,retmax,aliases[j],idtype)
					j+=1
			else: # add MESH term
				mesh=ref_mesh[i] # convert to original
				#mesh=mesh.split()
				mesh='+'.join(mesh.split()) # fill spaces    
				query=mesh+"[mesh] NOT review[Publication Type]"
				search[i]+=esearch_handle(db,retstart,retmax,query,idtype)
			# convert search to a series and append to searched 
			search=pd.Series(search,dtype=object)
			#print(search)
			searched=searched.append(search,verify_integrity=True)
			#print(searched)
	#
	## create combinations of "comb" keywords (default 4)
	def multinter(kws,searched=searched):
		ll=searched[list(kws)].tolist()
		#print(ll)
		ll=set(ll[0]).intersection(*ll)
		#print(ll)
		return ll 
	# at least 3 / at most "comb" entities must be used to formulate queries
	comb=min(comb,len(path))
	invc=len(path)-comb	
	hits=list(multinter(path[:comb])) # not imposing g1
	if invc != 0:
		hits+=list(multinter(path[invc:len(path)])) # not imposing g2
	hits=list(set(hits))
	#
	#print(hits)
	# retain only papers which are absent in the previous PMIDs list
	hits=[p for p in hits if p not in prev_ids]
	return hits,searched # spit out also the "searched" series
#
queryrecord=pd.DataFrame(columns=['path','pmids','len','cs'])
queryrecord.to_csv(wd+"gene-subgraphs_queries.tsv",sep='\t',mode='w',index=False)
#
def gapfind(s1,s2,comb=comb,gm=comp_Net,inter_query=inter_query,ref_mesh=ref_mesh,ref_gen=aliasdf):
	#
	global searched
	# s1 and s2 are dict items (key:net)
	# s1 will be a higher-priority network by construction
	k1,s1 = s1
	k2,s2 = s2
	#
	#print(s1.nodes(data=True))
	#
	s1=sorted(s1.nodes(data=True), key=lambda x: x[1]['ranking']) # list of (node,{attr})
	s2=sorted(s2.nodes(data=True), key=lambda x: x[1]['ranking']) # list of (node,{attr})
	#
	# we only need the node names
	s1=[t[0] for t in s1]
	s2=[t[0] for t in s2]
	#
	print("local removal external genes from network...")
	# ease computation time: delete gene nodes that belong to neither of the two subgraphs
	genes=[x for x,y in gm.nodes(data=True) if y['label']=="GENE"]
	rgenes=set(genes)-set(s1+s2)
	gm.remove_nodes_from(rgenes)
	#
	def pscore(g1,g2,s1=s1,s2=s2,Net=comp_Net):
		#
		ind1=s1.index(g1)+1
		ind2=s2.index(g2)+1
		if min([ind1,ind2])==ind1:
			# it has priority 
			prio=0
		else:
			# (1,2) > (2,1)
			prio=1
		#
		k=hmean([ind1,ind2])
		return ((g1,g2),k,prio)
	#
	# ranked list of gene pairs
	out=list(starmap(pscore,itertools.product(s1,s2)))
	out=sorted(out, key=itemgetter(1, 2))
	#out=out[:10]
	#
	npairs=iter([t[0] for  t in out])
	# now check for the first X node pairs s.t. > 100 papers can be retrieved, 
	# provided they contain the MeSH term connecting the genes
	# generate shortest paths
	def shrinkp(p,s1=s1,s2=s2):
		#print(p)
		# p is a list
		def GMseq(kw,s1=s1,s2=s2):
			if kw in s1:
				return 's1'
			elif kw in s2:
				return 's2'
			else:
				return 'M'
		#
		ls=list(map(GMseq,p))
		ls=''.join(ls) # a sequence like (regex) "s1+M+s2+"
		### ILLEGAL PATH: s1M...Ms1 ###
		if bool(re.search("s1Ms1|s2Ms2",ls)):
			return None
		#
		else:
			ls_out=re.split('s1[M]+s2',ls) #['',''] if no shrinkage is possible
			if len(ls_out)<2:
				ls_out=['','']
			#
			# if any s1/s2 is in-between, you would have smth like ['s1s1M','Ms2s2']
			#
			ls_out=[x.replace('M','') for x in ls_out]
			#
			#print(ls,ls_out)
			#
			s1_out= ls_out[0].count('s1') #len([i for i in ls_out if i=='s1'])
			s2_out=len(p) - ls_out[1].count('s2') #len([i for i in ls_out if i=='s2'])
			#
			#print(s1_out,s2_out)
			# e.g. ls_out==['s1','s2'] means exclude first gene and last gene of the path
			p=p[s1_out:s2_out]
			#print(p)
			return p
	#
	hpairs=[]
	#
	while (pair := next(npairs,False)):
		#
		g1,g2 = pair
		# search for all SHORTEST paths in the "complete" netw
		#
		print("get all shortest paths for %s-%s connections" % (g1,g2))
		l=list(nx.all_shortest_paths(gm,g1,g2,weight='inv_edge_weight'))
		#g1,g2 = pair
		# if a gene belonging to a graph is "closer", then restrict the query
		#	
		l=list(map(shrinkp,l))
		l=[i for i in l if i]
		#
		### do not allow queries which contain an inbetween gene not belonging to the two 
		aliases=set(aliasdf.symbol.to_list())
		okg=set(s1+s2)
		def checkg(g):
			if g in aliases and g not in okg: 
				return False # query must be discarded
			else:
				return True
		l1=[p for p in l if all([checkg(g) for g in p])]
		#
		#print(len(l1),len(l))
		#
		if any([len(p) > 0 and len(p) < 2*comb for p in l1]):
			hpairs.append(iter(l))
	#
	#############
	# limit the max number of paths to attempt a search
	#
	hpairs=hpairs[:P]
	#print(hpairs)
	#
	############	
	batch=[]
	retbatch=[]	
	for ps in hpairs:
		#
		ps=iter(ps)
		#print("try %s" % (str(ps)))
		#dataframe that records path-specific query outputs
		smbatch=pd.DataFrame(columns=['path','pmids','len'])	
		#
		while (path := next(ps,False)):
			#
			#
			print("search connections for %s" % (str(path)))
			#print('s\ntrying motifs pair:\n %s \n' % str(motp))
			#visited.append([str(mot) for mot in motp]) # append found motif to already visited motifs
			hits,searched = inter_query(path,searched,comb) # returns a tuple which is mapped to two objects
			hits=[h for h in hits if h not in batch]
			batch+=hits			
			batch=list(set(batch))
			#
			smbatch=smbatch.append(pd.DataFrame([['.'.join(path),hits,len(hits)]],columns=['path','pmids','len']))
			#
			# collect hits in a pandas series, and later sort from lowest to highest yield
			#
			print("%i paper(s) were found to join graphs %i and %i" % (len(batch),k1,k2))
			#
			#if >= 100 articles, we found enough for connecting graphs s1 and s2
			#if len(batch)>=100:
			#	break
			#
		#	
		smbatch.sort_values(by=['len','path'],inplace=True)
		# OC, THERE CAN STILL BE TIES
		smbatch.reset_index(drop=True,inplace=True)
		#
		smbatch['cs']=smbatch['len'].cumsum()
		#
		smbatch.to_csv(wd+"gene-subgraphs_queries.tsv",sep='\t',mode='a',index=False,header=False)
		#
		print("the queried paths yielded these lists of PMIDs:\n\n",smbatch,"\n\n")
		print("performing bottom-up cumulative gathering until all hits or a subset of at least %i PMIDs is reached" % (100))
		#
		if smbatch.len.sum() <= 100:
			# get all of them
			hits=smbatch.pmids.to_list()
			hits=list(set(itertools.chain.from_iterable(hits)))
			retbatch+=hits
		#
		else:
			# cumulative sum, then get the minimun index for a sum>=100
			#cs=smbatch['len'].cumsum().tolist()
			cs=smbatch.cs.to_list()
			trig=min([cs.index(p) for p in cs if p>=100])+1
			# drop everything above
			hits=smbatch.drop(list(range(trig,len(smbatch))))
			#print(hits)
			#
			hits=hits.pmids.to_list()
			hits=list(set(itertools.chain.from_iterable(hits)))
			retbatch+=hits
		#	
		print("\n### %i new paper(s) were selected to join graphs %i and %i ###\n" % (len(hits),k1,k2))
		#
		#if len(batch)>=100:
		#	break
	#					
	return retbatch			
#
#test=list(itertools.combinations(list(subgnets.items())[:3],2))[0]
#gapfind(*test)
#print(list(subgnets.items()))

print("start searching...")

batch=list(starmap(gapfind,itertools.combinations(list(subgnets.items()),2)))
batch=list(set(itertools.chain.from_iterable(batch)))
#
print("save searched queries")
compress_pickle(sd+'input/searched_queries',searched)
#
print("A total of %i new PMIDs were gathered to connect the subgraphs" % len(batch))

#
with open(prev_id_file,'a') as f:
	s=["\n\n", "######### SUB-GENE-GRAPHS JOINING ########",'\n']
	f.writelines(s)
	for pmid in batch:
		f.writelines('\n'+str(pmid))
#
if int(it) > 1:
	new_tag=re.sub("_subgraphs_expansion[\d]*", "_subgraphs_expansion"+it, tag)
	#new_tag=tag.split(r'_')[0]
	#print("new tag is (1): ",new_tag)
	#new_tag=new_tag+"_iteration"+it
	print("new tag is: ",new_tag)
else: # for the first iteration
	new_tag=tag+"_subgraphs_expansion1"#+it
#new_tag=tag+'' # specific tag
#
with open(tmp+"new_iteration.txt","w") as file:
	file.write(new_tag+"\n")

with open(tmp+"previous_iteration.txt","w") as file:
	file.write(tag+"\n")

# Save new pmids
with open(tmp+"efetch_inputs/"+new_tag+"_pmids_to_efetch.txt","w") as f:
	for pmid in batch:
		f.writelines(str(pmid)+'\n')

# save updated visited motifs
#visited=pd.DataFrame(visited)#,visited[0])
#visited.to_csv(tmp+"visited_motifs.tsv",sep='\t',index=False)
#
