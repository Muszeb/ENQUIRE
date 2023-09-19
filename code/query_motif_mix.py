# query mixtures of 2 motifs
# exec(open("/home/widmanmx/tam_textmining/code/query_motif.py").read())
# exec(open("code/query_motif_mix.py").read())
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
from scipy.stats import fisher_exact
from Bio import Entrez
from Bio.Entrez import efetch,esearch, read
from scipy.spatial import distance_matrix
import scipy
import csv
try:
	import _pickle as pickle
except ModuleNotFoundError:
	import pickle
#
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

#wd='/home/widmanmx/tam_textmining/'
#tag='iter0'
wd=sys.argv[1]
tag=sys.argv[2]
prev_id_file= sys.argv[3]
it=sys.argv[4]
comb=int(sys.argv[5])
sd=re.split("tmp",wd)[0]
tmp=re.split("\/"+tag,wd)[0]+'/'

# load motifs (gene1, gene2, mesh1, score) ranked highest to lowest score
#motif_tbl=pd.read_csv(wd+"data/"+tag+"_motifs_subs_and_ranks.tsv",delimiter='\t')
# QUERY QUEUE
qq=decompress_pickle(wd+"data/"+tag+"_query_queue.pbz2")
print("motifs loaded")
#
#aliasdf=pd.read_csv(sd+"input/gene_aliases_df.tsv",sep='\t')
aliasdf=pd.read_csv(sd+'input/aliases_uniprotensembl_as_df_Paralogues.tsv',sep='\t')
# remove bad aliases
aliases=aliasdf.alias_symbol.tolist()
# unless they are reference genes => regardless being a reference gene
#refs=aliasdf.symbol.tolist()
#
bad_hits=list(set(filter(lambda x: re.match("^[a-zA-Z]{1,2}$",str(x)),aliases)))
aliasdf=aliasdf[~aliasdf['alias_symbol'].isin(bad_hits)]
#defining functions
Entrez.email='Luca.musella@uk-erlangen.de'
Entrez.api_key='b89451c47164b8705c76108eec6386658607'
# GET #RETMAX PMIDS BY MESH TERM
# retstarts takes into account sorted list 
# retmax must be <= 100k
# transform mesh term back to original
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

# visited motifs (as pairs):
if int(it)<=1:
	visited=pd.read_csv(tmp+"visited_motifs.tsv",delimiter='\t',names=['Motif1','Motif2'])
else:
	visited=pd.read_csv(tmp+"visited_motifs.tsv",delimiter='\t')
visited=visited.values.tolist() # list of lists, the latter having string representation of tuples as values
# load the already searched queries, if any
try:
	searched=decompress_pickle(sd+'input/searched_queries.pbz2')
except:
	searched=pd.Series(dtype=object) # a series containing all KW (genes or MeSH) which have been extensively searched already
#
def gettype(t):
	if t in ref_gen.symbol.values:
		return t,'token'
	elif t in ref_mesh:
		return t,'mesh'
#	
def concmot(mm):
	return [list(set(t[0]+t[1])) for t in mm]  
#
Q={}
for q in qq: 
	Q.update(q)
#
#Q={k:concmot(v) for k,v in Q.items()}
#
# i=0
# # while i < len(terms):
# # 	Q=terms[i]
# # 	Q=[gettype(q) for q in Q]
# # 	res=list(starmap(query,Q))
# # #
# for i in range(0,len(terms)):
# 	Q=terms[i]
# 	Q=[gettype(q) for q in Q]
# 	res=list(starmap(query,Q))
# #	
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
print("formulate e-search queries")
#
def joinquery(g1g2):
	g1g2,qq = g1g2 # g1g2 is an element of Q.items()
	#
	print("attempt at joining graphs %s" % (str(g1g2)))
	#
	k=0
	batch=[] # new ID pool
	while k < len(qq):	# len qq being the number of attempts
		motp=qq[k]
		ents=[]
		if any([all([str(mot) in row for mot in motp]) for row in visited]):
			print('motifs pair already visited:\n %s \n' % str(motp))
			k+=1
		else:
			print('\ntrying motifs pair:\n %s \n' % str(motp))
			visited.append([str(mot) for mot in motp]) # append found motif to already visited motifs
			hits,searched = inter_query(motp[0]+motp[1]) # returns a tuple which is mapped to two objects
			batch+=hits
			batch=list(set(batch))
			# trigger: at least 100 new papers
			if len(hits)==0:
				print ("no new papers found for this motifs pair \n")
				k+=1
			elif len(batch)<100:			
				print('found',len(hits),' pubmed IDs for this motifs pair \n')
				print(len(batch),' new pubmed IDs have been gathered so far \n')
				ents+=list(motp[0]+motp[1])
				k+=1
			else: # trigger new iteration
			#write new id pool to previous id file
				with open(prev_id_file,'a') as f:
					s=["\n\n", "######### ITERATION "+it+"\t"+str(motp),'\n']
					f.writelines(s)
					for pmid in batch:
						f.writelines('\n'+str(pmid))
				# end big while loop
				print('\n',len(batch),'new pubmed IDs saved to file')
				Q=len(qq)
	return batch	
	#
### new papers to txt, update motif tsv
super_res={k:joinquery(i) for k,i in zip(Q.keys(),Q.items())}
#
batch=set(itertools.chain(*super_res.values()))
#
# reconstruct new tag
if int(it) > 1:
	new_tag=re.sub("_iteration[\d]+", "_iteration"+str(it), tag)
	#new_tag=tag.split(r'_')[0]
	#print("new tag is (1): ",new_tag)
	#new_tag=new_tag+"_iteration"+it
	print("new tag is: ",new_tag)
else: # for the first iteration
	new_tag=tag+"_iteration"+str(it)

# must be saved 1 level above the running tag-specific iteration
#tmp=re.split(tag+"\/",wd)[0]
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
visited=pd.DataFrame(visited)#,visited[0])
visited.to_csv(tmp+"visited_motifs.tsv",sep='\t',index=False)
#
# save searched queries
compress_pickle(sd+'input/searched_queries',searched)
# CREATE PRELIMINARY GENES AND MESH NETWORKS
sd=re.split("tmp",wd)[0]
tmp=re.split(tag+"\/"+"data",wd)[0]
exec(open(sd+"code/sort_network.py").read())