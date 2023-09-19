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
from os import listdir
from os.path import isfile, join
from collections import Counter
import copy
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
#
#wd='/home/widmanmx/tam_textmining/'
#tag="iter0"
#
wd= sys.argv[1]
# wd='/home/widmanmx/tam_textmining/'
tag= sys.argv[2] #iter1
#
thr=float(sys.argv[3])
#
sd=re.split("tmp",wd)[0]
tmp=''.join(re.split("\/"+tag,wd))
#
node_tbl=pd.read_csv(wd+tag+"_Genes_nodes_table_subgraph.tsv",delimiter='\t')
edge_tbl=pd.read_csv(wd+tag+"_Genes_edges_table_subgraph.tsv",delimiter='\t')
#
cols=list(node_tbl.columns)
# category/label is 2nd column
cat=cols[1]
#
attr={node_tbl[cols[0]][i]:{cols[j]:node_tbl[cols[j]][i] for j in range(1,len(cols))} for i in range(len(node_tbl))}
#
Network=nx.from_pandas_edgelist(edge_tbl,'kw1','kw2',['pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration'])
# adding node attributes
nx.set_node_attributes(Network,attr)
del attr
#
genes=node_tbl['KW'].tolist()
#
print("Network and attributes loaded")
# #
#dgenes={node_tbl.loc[i]['KW']:node_tbl.loc[i]['label'] for i in range(0,len(node_tbl))}
# def inmesh(g1,g2,G=Network,dg=dgenes):
# 	def t2g(l,dg):
# 		return(sum([dg[i]=='GENE' for i in l]))
# 	if not G.has_edge(g1,g2):
# 		try:
# 			asp=list(nx.all_shortest_paths(Network,g1,g2))
# 			asp=[sp for sp in asp if t2g(sp,dg)==2]
# 			return asp
# 		except nx.NetworkXNoPath:
# 			return []


# ### better parallelizing ###
# out=list(starmap(inmesh,itertools.combinations(genes,2)))
# #
# out=[o for o in out if o!= []] # each ELEMENT is still a list of lists
# structure of each element:
# [[g1,m1,m2,...,g2],[g1,mj,mk,...,g2]]
# DETECT GENE SUBGRAPHS AND STORE THEM SEPARATELY
print("detecting gene subgraphs...")
Net_genes=Network.subgraph(genes)
subs=list(Net_genes.subgraph(c).copy() for c in nx.connected_components(Net_genes))
# get graph size (nodes,edges)
# at least 3 nodes and 3 edges  -> product at least 9
subs={net:len(net.nodes())*len(net.edges()) for net in subs}
#
K=int(sys.argv[4])
#
nets=list(subs.keys())
for net in nets:
	if subs[net] < K:
		del subs[net]
#
ser_kwpub=decompress_pickle(wd+"data/"+tag+"_KWtoPub_series.pbz2") # IT HAS TO BE UPDATED ALONG THE ITERATION
# use m2h_map, if exists, to track forced conversion to human genes (absent in KWtoPub)
# try:
# 	m2h_map=decompress_pickle(tmp+"m2h_map_tie_nets.pbz2")
# 	#m2h_map=pd.Series(m2h_map)
# except:
# 	print("no mouse2human mapping found (not into sub-gene graphs expansion yet?)")
#

# 
tot_ids=len(set(itertools.chain.from_iterable(ser_kwpub.tolist())))
#
def represent(kw,tot_ids=tot_ids,ser=ser_kwpub):
	if kw in ser_kwpub.index:
		return len(ser[kw])/tot_ids
	else:
		print(kw)
		return 0
#
if len(subs)>=1:
	print("Found %i subgraph(s) with a 'characteristic connectivity' >= %i" % (len(subs),K))
	rank=list(subs.values())
	# create a sorted duplicate of rank, then translate to actual rankings
	srank=copy.deepcopy(rank)
	srank.sort(reverse=True)
	taken=[]
	for key in subs:
		r=srank.index(subs[key])#srank.index(subs[l[i]])#
		if r in taken:
			r=max(taken)+1	
		taken.append(r)
		subs[key]=r
	# create separated edge and table list
	print('calculating centralities (sub-net specific)')
	for key,val in subs.items():
		#
		# PRUNE BY REPRESENTATIVENESS
		nodes=list(key.nodes())
		# convert if a m2h map is available
		# if "m2h_map" in globals():
		# 	nodes=m2h_map[nodes].tolist()
		# #
		repGT=100*(len(set(itertools.chain.from_iterable(ser_kwpub[ser_kwpub.index.intersection(nodes)].tolist())))/tot_ids)
		repGM=100*(hmean(list(map(represent,nodes))))
		print('graph:',val,repGT, repGM)
		if repGT >= thr: # pruning threshold
			#
			# closeness=nx.closeness_centrality(key,distance='inv_edge_weight')
			# nx.set_node_attributes(key,closeness,"closeness_cnt")
			# #
			# betweenness=nx.betweenness_centrality(key,weight='inv_edge_weight')
			# nx.set_node_attributes(key,betweenness,"betweenness_cnt")
			# #
			# strength=dict(key.degree(weight='edge_weight'))
			# nx.set_node_attributes(key,strength,'strength')
			#
			edge_tbl=nx.to_pandas_edgelist(key,source='kw1',target='kw2')
			node_tbl=pd.DataFrame.from_dict(dict(key.nodes(data=True)), orient='index')
			node_tbl=node_tbl.reset_index()
			node_tbl=node_tbl.rename(columns={"index": "entity"})
			node_tbl['repGM']=[repGM]*len(node_tbl)
			node_tbl['repGT']=[repGT]*len(node_tbl)
			node_tbl.to_csv(wd+"/data/gene-subgraphs/"+tag+"_nodes_table_gene-subgraph_"+str(val)+".tsv",sep='\t',index=False)
			edge_tbl[['kw1','kw2','pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration']].to_csv(wd+"/data/gene-subgraphs/"+tag+"_edges_table_gene-subgraph_"+str(val)+".tsv",sep='\t',index=False)

		#
print("gene-only subgraphs stored")