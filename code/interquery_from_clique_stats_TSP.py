import os
import sys
import numpy as np 
import pandas as pd
import itertools
from itertools import islice,product,starmap
import multiprocessing as mp
from collections import OrderedDict, Counter
import re
import random
from random import sample
import time
import networkx as nx
import bz2
import json
import statistics
from multiprocess import Process, Manager
from scipy.stats import fisher_exact
from Bio import Entrez
from Bio.Entrez import efetch,esearch, read
from scipy.spatial import distance_matrix
import networkx.algorithms.approximation as nx_app
import math
import scipy
import csv
from itertools import islice
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
print("load metadata...")
#wd='/home/widmanmx/tam_textmining/'
#
wd=sys.argv[1]
tag=sys.argv[2]
prev_id_file= sys.argv[3]
query_file=sys.argv[5]
webenvfile=tmp+tag+"_EDirect_WebEnv.txt"
#
sd=re.split("tmp",wd)[0]
tmp=re.split("\/"+tag,wd)[0]+'/'
#
it=sys.argv[4]
comb=int(sys.argv[6])
A=int(sys.argv[7])
#
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

# load the already searched queries, if any
try:
	searched=decompress_pickle(sd+'input/searched_queries.pbz2')
except:
	searched=pd.Series(dtype=object) # a series containing all KW (genes or MeSH) which have been extensively searched already
#
####
print("##############")
print("Start esearch queries...")
####
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
####
handle = Entrez.esearch(db="pubmed", term="lung+cancer")
record = Entrez.read(handle)
count = record['Count']
handle = Entrez.esearch(db="pubmed", term="lung+cancer", retmax=count)
record = Entrez.read(handle)
request.post("esearch.fcgi?db=<database>&term=protein>")
# catch papers above the retmax limit
while count>retstart:
	try:
		print("search for %s: " % (term),count, retstart)
		#print("search for %s",count,retstart) % (term)
		handle = Entrez.esearch(db=db,retstart=retstart, retmax=retmax, term=term, idtype=idtype)
		record = Entrez.read(handle)
		handle.close()
		res += record['IdList']
		count=int(record['Count']) # actually immutable
		
	except Exception as e:
		print(e)
		print("Internet connection problems?")
		time.sleep(5)
	#
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
			if  i.count('/') == 1:
				mesh,subh= i.split('/') # take care of subheadings
				mesh=ref_mesh[mesh] # convert to original
				#mesh=mesh.split()
				mesh='+'.join(mesh.split()) # fill spaces
				subh='+'.join(subh.split()) # fill spaces
				query=mesh+'/'+subh+"[mesh] NOT review[Publication Type]"
			elif i.count('/') == 0:
				mesh=ref_mesh[i] # convert to original
				mesh='+'.join(mesh.split()) # fill spaces
				query=mesh+"[mesh] NOT review[Publication Type]"
			else: 
				print("WARNING! MULTIPLE SUBHEADINGS DETECTED! CHECK CODE")
				mesh,subh= i.split('/')[:2] # take care of subheadings
				mesh=ref_mesh[mesh] # convert to original
				#mesh=mesh.split()
				mesh='+'.join(mesh.split()) # fill spaces
				subh='+'.join(subh.split()) # fill spaces
				query=mesh+'/'+subh+"[mesh] NOT review[Publication Type]"
			#
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
###
conqnet=pd.read_csv(query_file,sep='\t')
### edges should be the smallest when the connection is strongest
conqnet['weight']=conqnet.weight.apply(lambda x: 1/(x+0.1)) # handle 0s 
### TRAVELLING SALESMAN PROBLEM ###
tsp = nx.approximation.traveling_salesman_problem
# will use christofides as default #
### check how pairs repeat themselves ### 
def window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result
# ###
# pairs=[]
# ##
# for i in range(10000):
# 	conqnet=conqnet.sample(frac=1)
# 	best_p=tsp(net)
# 	pairs.append(list(set(window(best_p))))
# ###
# pairs=list(itertools.chain.from_iterable(pairs))
# # >>> pairs=list(itertools.chain.from_iterable(pairs))
# # >>> from collections import OrderedDict, Counter
# # >>> Counter(pairs)
# # Counter({(7, 4): 10000, (8, 4): 10000, (2, 7): 10000, (3, 1): 10000, (1, 8): 10000, (4, 6): 10000, (6, 4): 10000, (4, 2): 10000, (4, 5): 10000, (5, 3): 10000})
# # >>> 
# ### VERY STABLE! ###
#

### we will have various counters: ###
#
# 1) for PMIDs (batch)
# 2) for the visited queries (visited)
# 3) for the communities for which a query has successfully retrieved 
#
####
batch=[]
visited=[]
commconn=[]
communconn=[]
#
j=0
while j < A: 
	#
	print("### ATTEMPT N. %i ###" %(j+1))
	#
	net=nx.from_pandas_edgelist(conqnet[~conqnet['query'].isin(visited)],'c1','c2',['weight','query'])
	#### best TSP ####
	best_p=tsp(net)
	print("approximately best TSP path is %s \n" % (str(best_p)))
	### get sliding window and single pairs ###
	best_pairs=set([tuple(sorted(list(p))) for p in window(best_p)])
	# assuming the smallest of multi-edge options were used... #
	def candiq(t,conqnet=conqnet):
		#print(t)
		c1,c2=t
		# idxmin min guarantees to return one single value #
		bestq=conqnet.loc[conqnet[((conqnet['c1']==c1) & (conqnet['c2']==c2)) | ((conqnet['c1']==c2) & (conqnet['c2']==c1))].weight.idxmin()].query
		return bestq
	#
	best_qs={p:candiq(p) for p in best_pairs}
	### add to visited ###
	visited=visited+list(best_qs.values())
	###
	best_qs={k:v.split('_') for k,v in best_qs.items()}
	print("Best Candidate queries:")
	print(best_qs)
	print("\n")
	#### BLAME PUBMED UPDATE OF ESEARCH FROM 21.11.22 ####
	print("Annotate query keys and generate WebEnv via EDirect ...")
	kws=frozenset(itertools.chain.from_iterable(best_qs.values()))
	querykeys=pd.Series({k:"(#%s)" % (i) for k,i in zip(kws,range(1,len(kws)+1))})
	#
	### format based on mesh/gene ###
	def formatkw(i):
		if i in aliasdf.symbol.values: # it's a gene							
			# detect then all aliases
			aliases=aliasdf[aliasdf['symbol']==i]['alias_symbol'].tolist()
			aliases=['+'.join(al.split())+"[Title/Abstract]" for al in aliases] # add "NOT review" at the end NOT review[Publication Type]" for al in aliases]
			return " OR ".join(aliases)
		else: # it's a MESH term
			if  i.count('/') == 1:
				mesh,subh= i.split('/') # take care of subheadings
				mesh=ref_mesh[mesh] # convert to original
				#mesh=mesh.split()
				mesh='+'.join(mesh.split()) # fill spaces
				subh='+'.join(subh.split()) # fill spaces
				query=mesh+'/'+subh+"[Mesh]" # NOT review[Publication Type]"
			elif i.count('/') == 0:
				mesh=ref_mesh[i] # convert to original
				mesh='+'.join(mesh.split()) # fill spaces
				query=mesh+"[Mesh]" #" NOT review[Publication Type]"
			else: 
				print("WARNING! MULTIPLE SUBHEADINGS DETECTED! CHECK CODE")
				mesh,subh= i.split('/')[:2] # take care of subheadings
				mesh=ref_mesh[mesh] # convert to original
				#mesh=mesh.split()
				mesh='+'.join(mesh.split()) # fill spaces
				subh='+'.join(subh.split()) # fill spaces
				query=mesh+'/'+subh+"[Mesh]" # NOT review[Publication Type]"
			return query

	kws2q=frozenset({formatkw(i) for i in kws})
	kws2q=frozenset({"'" + i + "'" for i in kws})
	### additionally, construct already the necessary combinations ###
	interqueries={k:" AND ".join(querykeys[v].tolist()) for k,v in best_qs.items()}
	interqueries={k:"'("+v+") NOT (review[PublicationType])'" for k,v, in interqueries.items()}
	interkeys=pd.Series({k:"(#%s)" % (i) for k,i in zip(best_qs.keys(),range(len(querykeys)+1,len(querykeys)+1+len(best_qs)+1))})
	### construct piped esearch and save web environment ### 
	pipe="esearch -db pubmed -query "+" | esearch -db pubmed -query ".join(kws2q)
	pipe2=" | esearch -db pubmed -query "+" | esearch -db pubmed -query ".join(interqueries.values())
	###
	os.system("cat "+pipe+pipe2+" > "+wd+tag+"_EDirect_WebEnv.txt")
	####


#
for k,v in best_qs.items():
	if any([com not in commconn for com in k]):
		#
		print("try connecting communities %s and %s ..." % (k[0],k[1]))
		print("search for  %s" % (str(v)))
		hits,searched = inter_query(v,searched,comb)
		#
		if len(hits)>0:
			commconn=list(set(commconn+list(k)))
			communconn=list(set(communconn)-set(commconn))
		else: 
			communconn=list(set(communconn+list(k)))
		#
		batch=list(set(batch+hits))
		print("A total of %i new PMIDs has been collected \n" % (len(batch)))
#
	if len(communconn)==0:
		print("All communities produced a successful query \n")
		break
	else:
		print("communities %s could not produce a successful query \n" % (str(communconn)))
		j+=1	
	#

####
resQ=batch
####
print("save searched queries")
compress_pickle(sd+'input/searched_queries',searched)
#
with open(prev_id_file,'a') as f:
	s=["\n\n", "######### SUB-GENE-GRAPHS JOINING ########",'\n']
	f.writelines(s)
	for pmid in resQ:
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
	for pmid in resQ:
		f.writelines(str(pmid)+'\n')
#
