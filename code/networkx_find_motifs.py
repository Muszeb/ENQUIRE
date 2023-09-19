# networkx test
# exec(open("/home/widmanmx/tam_textmining/code/networkx_find_motifs.py").read())
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
from collections import Counter
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

### import Network after ranking nodes in R ###

edge_tbl=pd.read_csv(wd+"data/"+tag+"_allxall_sig_combinations_bh_to_cytoscape.tsv",delimiter='\t')
ew_list=list(edge_tbl['edge_weight'])
edge_tbl['1-edge_weight']=[1-ew for ew in edge_tbl['edge_weight']]

node_tbl=pd.read_csv(wd+"data/"+tag+"_node_table_with_rankings_strength_only.tsv",delimiter='\t')


#node_tbl=pd.read_csv(wd+"data/"+tag+"_node_table_with_rankings.tsv",delimiter='\t')


#converting node_tbl to arrtibutes dict
attr={node_tbl['label'][i]:{
'category':node_tbl['category'][i],
'weight':node_tbl['node_weight'][i],
'cls_cnt':node_tbl['cls_cnt'][i],
'btw_cnt':node_tbl['btw_cnt'][i],
'strength':node_tbl['strength'][i],
'Krank':node_tbl['Krank'][i],
'ranking':node_tbl['ranking'][i]}
for i in range(len(node_tbl))}

# creating Network
Network=nx.from_pandas_edgelist(edge_tbl,'kw1','kw2',['edge_weight','1-edge_weight'])

# adding node attributes
nx.set_node_attributes(Network,attr)
del attr
#
print('Network loaded')
#
### finding motif GENE-GENE-MESH closed triplet ###

def get_genelist(node,g):
	l=[]
	if node in g:
		for nb in g.neighbors(node):
			if g.nodes('category')[nb]=='GENE':
				l.append(nb)

		if len(l)>=2:
			gene_nb=(node,l)
			return gene_nb

# save as tuples gene_nb = (origin mesh,[list of genes])


#
def get_nodes(n1,n2,g=Network):
	if g.has_edge(n1,n2):
		return (n1,n2)

# to return nodes instead of bool with starmap

#
def find_motif(mesh_gene_list_tuple,all_motifs):
	motifs=[]
	genes=mesh_gene_list_tuple[1]
	edge_list=[x for x in starmap(get_nodes,itertools.combinations(genes,2)) if x]
	if edge_list:
		for e in edge_list:
			motifs.append((e[0],e[1],mesh_gene_list_tuple[0]))
		if motifs:
			all_motifs.append(motifs)

# returns motif as list of tuples(gene1,gene2,origin mesh)

mesh_list=[node for node in list(Network.nodes) if Network.nodes('category')[node]=='MESH']

all_motifs=[]

for mesh in mesh_list:
	g_genes=(get_genelist(mesh,Network))
	if g_genes:
		find_motif(g_genes,all_motifs)

all_motifs=list(itertools.chain.from_iterable(all_motifs))

all_motifs=pd.DataFrame(all_motifs, columns=['gene1','gene2','mesh'])

# rank by K-node score 

# def add_scores(i,df=all_motifs):
# 	score=sum([Network.nodes('Krank')[node] for node in df.loc[i]])
# 	return score

# all_motifs['score']=[add_scores(i) for i in range(len(all_motifs))]

# all_motifs.sort_values(by='score',ascending=False,inplace=True)
# all_motifs.reset_index(drop=True,inplace=True)

# rank by ranking

def harm_mean(i, df=all_motifs):
	ranks=[Network.nodes('ranking')[node] for node in df.loc[i]]
	score=hmean(ranks)
	return score


all_motifs['score']=[harm_mean(i) for i in range(len(all_motifs))]

all_motifs.sort_values(by='score',ascending=True,inplace=True)
all_motifs.reset_index(drop=True,inplace=True)



all_motifs.to_csv(wd+"data/"+tag+"_motifs_strengthonly.tsv",sep='\t')

print('motifs saved')