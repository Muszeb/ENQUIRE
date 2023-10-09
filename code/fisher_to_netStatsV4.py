# CREATE NETWORK FROM DICT
# exec(open("/home/musellla/tam_textmining/code/create_network_from_dict.py").read())
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
from scipy.stats import rankdata
import threading
from time import sleep
from datetime import datetime
import scipy
import csv
from datetime import datetime
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
wd=sys.argv[1] #'/home/musellla/tam_textmining/'
tag=sys.argv[2] #"iter1"
thr=float(sys.argv[4])
#
edge_tbl=pd.read_csv(wd+"data/"+tag+"_allxall_sig_combinations_bh_to_cytoscape.tsv",delimiter='\t')
node_tbl=pd.read_csv(wd+"data/"+tag+"_allxall_pre_nodes_table_bh_to_cytoscape.tsv",delimiter='\t')
#
node_tbl=node_tbl.drop_duplicates(subset=['KW'])
edge_tbl=edge_tbl.drop_duplicates(subset=['kw1','kw2'])
# you have to consider swapped occurrences too
oks=pd.DataFrame()
oks['kw1']=edge_tbl[['kw1','kw2']].apply(lambda x: x.sort_values().iloc[0],axis=1)
oks['kw2']=edge_tbl[['kw1','kw2']].apply(lambda x: x.sort_values().iloc[1],axis=1)
oks=list(oks.drop_duplicates().index)
edge_tbl=edge_tbl.filter(oks,axis=0)
#
edge_tbl.reset_index(drop=True,inplace=True)
node_tbl.reset_index(drop=True,inplace=True)
#
edge_tbl.to_csv(wd+"data/"+tag+"_allxall_sig_combinations_bh_to_cytoscape.tsv",sep='\t',index=False)
#
print("compute Jaccard distance-weighted network parameters: \n 1) closeness centrality \n 2) betweenness centrality \n 3) strength")

l_dict={node_tbl['KW'][i]:node_tbl['category'][i] for i in range(len(node_tbl))}

w_dict={node_tbl['KW'][i]:node_tbl['node_weight'][i] for i in range(len(node_tbl))}

o_dict={node_tbl['KW'][i]:node_tbl['occurrence'][i] for i in range(len(node_tbl))}
#
#network=nx.from_pandas_edgelist(edge_tbl,'kw1','kw2',['OR','pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration'])
network=nx.from_pandas_edgelist(edge_tbl,'kw1','kw2',['pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration'])
#
nx.set_node_attributes(network,l_dict,"category")
nx.set_node_attributes(network,w_dict,"node_weight")
nx.set_node_attributes(network,o_dict,"occurrence")
#
print('Network loaded')
#
print('calculating centralities (whole network)')
#
# closeness=nx.closeness_centrality(network,distance='inv_edge_weight')
# nx.set_node_attributes(network,closeness,"closeness_cnt")
# #
# betweenness=nx.betweenness_centrality(network,weight='inv_edge_weight')
# nx.set_node_attributes(network,betweenness,"betweenness_cnt")
# #
# #strength=dict(network.degree(weight='edge_weight'))
# strength=dict(network.degree(weight='R1C1'))
# nx.set_node_attributes(network,strength,'strength')
# #
# DETECT SUBGRAPHS AND STORE THEM SEPARATELY
print("detecting subgraphs...")
subs=list(network.subgraph(c).copy() for c in nx.connected_components(network))
# get graph size (nodes,edges)
# at least 3 nodes and 3 edges  -> product at least 9
subs={net:len(net.nodes())*len(net.edges()) for net in subs}
#
nets=list(subs.keys())
for net in nets:
	if subs[net] < int(sys.argv[3]): 
		del subs[net]
#
if len(subs)>=1:
	print("Found %i subgraphs with >= 3 nodes, >= 3 edges" % (len(subs)))
	rank=list(subs.values())
	srank=copy.deepcopy(rank)
	srank.sort(reverse=True)
	taken=[]
	#
	ser_kwpub=decompress_pickle(wd+"data/"+tag+"_KWtoPub_series.pbz2") # IT HAS TO BE UPDATED ALONG THE ITERATION
	tot_ids=len(set(itertools.chain.from_iterable(ser_kwpub.tolist())))
	#
	def represent(kw,tot_ids=tot_ids,ser=ser_kwpub):
		return len(ser[kw])/tot_ids
	#
	print("compute representation parameters...")
	#
	for key in subs:
		r=srank.index(subs[key])#srank.index(subs[l[i]])#
		if r in taken:
			r=max(taken)+1	
		taken.append(r)
		subs[key]=r
	# create separated edge and table list
	k=0
	for key,val in subs.items():
		#print('calculating centralities (sub-net specific)')
		#
		# closeness=nx.closeness_centrality(key,distance='inv_edge_weight')
		# nx.set_node_attributes(key,closeness,"closeness_cnt")
		# #
		# betweenness=nx.betweenness_centrality(key,weight='inv_edge_weight')
		# nx.set_node_attributes(key,betweenness,"betweenness_cnt")
		# #
		# strength=dict(key.degree(weight='R1C1'))
		# nx.set_node_attributes(key,strength,'strength')
		#
		edge_tbl=nx.to_pandas_edgelist(key,source='kw1',target='kw2')
		node_tbl=pd.DataFrame.from_dict(dict(key.nodes(data=True)), orient='index')
		node_tbl=node_tbl.reset_index()
		node_tbl=node_tbl.rename(columns={"index": "entity"})
		repGT=100*(len(set(itertools.chain.from_iterable(ser_kwpub[node_tbl.entity.tolist()].tolist())))/tot_ids)
		if repGT >= thr:		
			k+=1
			node_tbl.to_csv(wd+"data/subgraphs/"+tag+"_nodes_table_subgraph_"+str(val)+".tsv",sep='\t',index=False)
			edge_tbl.to_csv(wd+"data/subgraphs/"+tag+"_edges_table_subgraph_"+str(val)+".tsv",sep='\t',index=False)
		#
		print("%i subgraph(s) passed the %g%% threshold" % (k,thr))
	#

#
#save Network
print("save nodes table")
new_node_tbl=pd.DataFrame.from_dict(dict(network.nodes(data=True)), orient='index')
new_node_tbl=new_node_tbl.reset_index()
new_node_tbl=new_node_tbl.rename(columns={"index": "KW"})
#
new_node_tbl.to_csv(wd+"data/"+tag+"_allxall_nodes_table_bh_to_cytoscape.tsv",sep='\t',index=False)