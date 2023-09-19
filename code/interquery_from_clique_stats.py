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
import statistics
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
qq=pd.read_csv(query_file,sep='\t')
#
qq=qq.apply(lambda x : [v.split('_') for v in x.tolist() if type(v)==str],axis=1).to_dict()
#
print("...done")
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
				if  i.count('/') == 1:
					mesh,subh= i.split('/') # take care of subheadings
					mesh=ref_mesh[mesh] # convert to original
					#mesh=mesh.split()
					mesh='+'.join(mesh.split()) # fill spaces
					subh='+'.join(subh.split()) # fill spaces
					query=mesh+'/'+subh"[mesh] NOT review[Publication Type]"
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
					query=mesh+'/'+subh"[mesh] NOT review[Publication Type]"
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
#
def triggerQ(k):
	global searched
	batch=[]
	for q in qq[k]:	
		print("search for  %s" % (str(q)))
		hits,searched = inter_query(q,searched,comb)
		
		batch=list(set(batch+hits))
		print("A total of %i new PMIDs has been collected" % (len(batch)))
		if len(batch) != 0:
			break
	return batch
#
resQ=[triggerQ(k) for k in qq]
resQ=list(set(itertools.chain.from_iterable(resQ)))
#
print("save searched queries")
compress_pickle(sd+'input/searched_queries',searched)
#
print("A total of %i new PMIDs were gathered to connect the subgraphs" % len(resQ))

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
