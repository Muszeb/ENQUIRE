import os
import sys
import numpy as np 
import pandas as pd
import subprocess
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
prev_ids=[s.split('\n')[0] for s in prev_ids]
prev_ids=list(set(prev_ids))

# # load the already searched queries, if any
# try:
# 	searched=decompress_pickle(sd+'input/searched_queries.pbz2')
# except:
# 	searched=pd.Series(dtype=object) # a series containing all KW (genes or MeSH) which have been extensively searched already
# #
# ####
print("##############")
####
conqnet=pd.read_csv(query_file,sep='\t')
### edges should be the smallest when the connection is strongest
conqnet['weight']=conqnet.weight.apply(lambda x: x+0.1) # handle 0s 
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

##### QUERY EXAMPLE IN EDIRECT ##### 
# esearch -db pubmed -query 'anti-bacterial agents[mesh]' | \
# esearch -db pubmed -query 'clostridium infections/chemically induced[mesh]' | \
# esearch -db pubmed -query 'diarrhea[mesh]' | \
# esearch -db pubmed -query 'feces[mesh]' | \
# esearch -db pubmed -query '(#1) AND (#2) AND (#3) AND (#4) NOT review[Publication Type]'
### we will have various counters: ###
#
# 1) for PMIDs (batch)
# 2) for the visited queries (visited)
# 3) for the communities for which a query has successfully retrieved 
#
#
#### FUNCTIONS #####
def formatkw(i):
	if i in aliasdf.symbol.values: # it's a gene							
		# detect then all aliases
		aliases=aliasdf[aliasdf['symbol']==i]['alias_symbol'].tolist()
		aliases=['+'.join(al.split())+"[Title/Abstract]" for al in aliases] # add "NOT review" at the end NOT review[Publication Type]" for al in aliases]
		return "("+" OR ".join(aliases)+")"
	else: # it's a MESH term
		if  i.count('/') == 1:
			mesh,subh= i.split('/') # take care of subheadings
			mesh=ref_mesh[mesh] # convert to original
			#mesh=mesh.split()
			mesh='+'.join(mesh.split()) # fill spaces
			subh='+'.join(subh.split()) # fill spaces
			query=mesh+'/'+subh+"[mesh]" # NOT review[Publication Type]"
		elif i.count('/') == 0:
			mesh=ref_mesh[i] # convert to original
			mesh='+'.join(mesh.split()) # fill spaces
			query=mesh+"[mesh]" #" NOT review[Publication Type]"
		else: 
			print("WARNING! MULTIPLE SUBHEADINGS DETECTED! CHECK CODE")
			mesh,subh= i.split('/')[:2] # take care of subheadings
			mesh=ref_mesh[mesh] # convert to original
			#mesh=mesh.split()
			mesh='+'.join(mesh.split()) # fill spaces
			subh='+'.join(subh.split()) # fill spaces
			query=mesh+'/'+subh+"[mesh]" # NOT review[Publication Type]"
		return query

###
def uidlist(k,keys):
	keysform=kws2q[keys].values.tolist()
	#keysform=["'"+v+"'" for v in keysform]
	print("try connecting communities %s and %s ...\n" % (k[0],k[1]))
	# efetch 1st kw, efilter with the others, then filter out reviews #
	query="esearch -db pubmed -query '" + " AND ".join(keysform)+" NOT review[Publication Type]'"#" efilter -query " % (webenvfile,keysform[0])
	#query="cat %s | esearch -query %s | efilter -query " % (webenvfile,keysform[0])
	#query1=" | efilter -query ".join(keysform[1:])
	query=query+" | efetch -format uid"
	print(query)
	# shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT
	hits=subprocess.Popen(query,shell=True,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)#.read().split('\n')
	hits=hits.communicate()[0].decode('utf-8').split('\n')
	#hits=os.popen("cat test_out.txt | esearch -query '%s' | efetch -format uid" % (key)).read().split('\n')
	hits=list(set(hits)-set(prev_ids))
	print("\nThe query generated %i new PMIDs\n" % (len(hits)-1))
	return hits

####
batch=[]
visited=[]
commies=set(conqnet.c1.tolist()+conqnet.c2.tolist())
commconn=set()
communconn=set()
####
j=0
queryrecord=dict()
while j < A: 
	#
	print("### ATTEMPT N. %i ###" %(j+1))
	#
	net=nx.from_pandas_edgelist(conqnet[~conqnet['query'].isin(visited)],'c1','c2',['weight','query'])
	###
	### The Communities might belong to different subgraphs! ###
	#### best TSP ####
	S = [net.subgraph(c).copy() for c in nx.connected_components(net)]
	best_p=[tsp(n) for n in S]
	#if len(best_p)==1:
	#	best_p=[best_p]
	#
	#best_p=tsp(net)
	print("approximately best TSP paths are %s \n" % (' AND '.join([str(p) for p in best_p])))
	### get sliding window and single pairs ###
	best_pairs=[set([tuple(sorted(list(p))) for p in window(paths)]) for paths in best_p]
	best_pairs=set.union(*best_pairs)
	### filter based pairs based on still-unconnected communities ### 
	best_pairs={t for t in best_pairs if any([x not in commconn for x in t])}
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
	#best_qs={(4,7):best_qs[(4,7)]}
	print("Best Candidate queries:")
	print(best_qs)
	print("\n")
	#### BLAME PUBMED UPDATE OF ESEARCH FROM 21.11.22 ####
	print("Format query keys...")
	kws=frozenset(itertools.chain.from_iterable(best_qs.values()))
	querykeys=pd.Series({k:"(#%s)" % (i) for k,i in zip(kws,range(1,len(kws)+1))})
	#
	### format based on mesh/gene ###
	#
	#kws2q=frozenset({formatkw(i) for i in kws})
	kws2q=pd.Series({k:formatkw(k) for k in kws})
	#kws2q=frozenset({"'" + i + "'" for i in kws2q})
	#kws2form={k:f for k,f in zip(kws,kws2q)}
	### additionally, construct already the necessary combinations ###
	#interqueries={k:" AND ".join(querykeys[v].tolist()) for k,v in best_qs.items()}
	# 
	#interqueries={k:"'"+v+"'" for k,v in interqueries.items()}
	#interqueries={k:"'"+v+" NOT review[Publication Type]'" for k,v, in interqueries.items()}
	#interkeys=pd.Series({k:"(#%s)" % (i) for k,i in zip(best_qs.keys(),range(len(querykeys)+1,len(querykeys)+1+len(best_qs)+1))})
	### construct piped esearch and save web environment ### 
	#pipe="| esearch -db pubmed -query "+" | esearch -db pubmed -query ".join(kws2q)
	#pipe2=" | esearch -db pubmed -query "+" | esearch -db pubmed -query ".join(interqueries.values())
	###
	#os.system("cat " + webenvfile + pipe + " > " + webenvfile)
	####s
	#print("WebEnv generated, extract PMID lists using efetch...")
	#
	#
	fetchpmids={k:uidlist(k,v) for k,v in best_qs.items()}
	fetchpmids={k:[i for i in v if i != ''] for k,v in fetchpmids.items()}
	queryrecord.update({tuple(best_qs[k]):v for k,v in fetchpmids.items()})	
	#
	if all([len(i)==0 for i in fetchpmids.values()]):
		print("WARNING: all queries from ATTEMPT N.%i were unsuccessful" % (j+1))
		j+=1
	elif all([len(i)!=0 for i in fetchpmids.values()]):
		print("All communities produced a successful query \n")
		hits=list(set(itertools.chain.from_iterable(fetchpmids.values())))
		batch=list(set(batch+hits))
		print("A total of %i new PMIDs has been collected \n" % (len(batch)))	
		break
	else: 
		commconn=commconn.union(set(itertools.chain.from_iterable([k for k in fetchpmids.keys() if len(fetchpmids[k])>0])))	
		communconn=set(commies)-set(commconn)
		print("Communities %s did not produce a successful query\n" % (str(communconn)))
		hits=list(set(itertools.chain.from_iterable(fetchpmids.values())))
		batch=list(set(batch+hits))
		print("A total of %i new PMIDs has been collected \n" % (len(batch)))
		j+=1
# #
# for k,v in best_qs.items():
# 	if any([com not in commconn for com in k]):
# 		#
# 		print("try connecting communities %s and %s ..." % (k[0],k[1]))
# 		print("search for  %s" % (str(v)))
# 		hits,searched = inter_query(v,searched,comb)
# 		#
# 		if len(hits)>0:
# 			commconn=list(set(commconn+list(k)))
# 			communconn=list(set(communconn)-set(commconn))
# 		else: 
# 			communconn=list(set(communconn+list(k)))
# 		#
# 		batch=list(set(batch+hits))
# 		print("A total of %i new PMIDs has been collected \n" % (len(batch)))
# #
# 	if len(communconn)==0:
# 		print("All communities produced a successful query \n")
# 		break
# 	else:
# 		print("communities %s could not produce a successful query \n" % (str(communconn)))
# 		j+=1	
# 	#

####
resQ=batch
####
print("save searched queries")
#compress_pickle(sd+'input/searched_queries',searched)
#
with open(prev_id_file,'a') as f:
	s=["\n\n", "######### COMMUNITY-JOINING N.%s ########" % (it),'\n']
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
# Save query-to-PMID information #
pmids=list(itertools.chain.from_iterable(queryrecord.values()))
queries=list(itertools.chain.from_iterable([[k]*len(v) for k,v in queryrecord.items()]))
queryrecord={i:j for i,j in zip(pmids,queries)}
queryrecdf=pd.DataFrame(queryrecord.items(), columns=['PMID','QuerySet'])
queryrecdf.to_csv(tmp+"efetch_inputs/QueryToPMID_record.tsv",sep='\t',index=False,mode='a', header=False)