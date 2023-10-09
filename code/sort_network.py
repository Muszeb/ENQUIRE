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
# if 'sys.argv' in globals():
# 	if len(sys.argv)>1:
#
try:
	#print(sys.argv)
	wd=sys.argv[1]
	tag=sys.argv[2]
	sd=re.split("tmp",wd)[0]
	tmp=re.split(tag+"\/"+"data",wd)[0]
except: 
	print("no additional command-line arguments found")		
#	
#wd='/home/widmanmx/tam_textmining/'
#tag="iter0"
#
#thr= int(sys.argv[3])
#
edge_tbl=pd.read_csv(wd+"data/"+tag+"_allxall_sig_combinations_bh_to_cytoscape.tsv",delimiter='\t')#ew_list=list(edge_tbl['edge_weight'])
#edge_tbl['1-edge_weight']=[1-ew for ew in edge_tbl['edge_weight']]
### GENERATE HYPERLINKS TO PAPERS ###
pubmed='https://pubmed.ncbi.nlm.nih.gov/'
hd=re.split("\/"+tag,wd)[0]+'/'
kw2pubs=decompress_pickle(wd+"data/"+tag+"_KWtoPub_series.pbz2")
### 
def hyperlink(v):
	l=kw2pubs[v].tolist()
	l=list(set(l[0]).intersection(*l))
	l=[str(i) for i in l]
	def link(s):
		return '=HYPERLINK("'+pubmed+s+'/", "'+s+'")'
	links=[link(s) for s in l]
	links_expl=[pubmed+s for s in l]
	return [{'PMID': i, 'kw1': v.iloc[0], 'kw2': v.iloc[1], 'link':s, 'link_explicit':t} for i,s,t in zip(l,links,links_expl)]
###	
hypes=edge_tbl[['kw1','kw2']].apply(hyperlink,axis=1).tolist()
hypes=list(itertools.chain.from_iterable(hypes))
hypes=pd.DataFrame(hypes)
#
### JOIN WITH QUERY ATTRIBUTE ### 
tmp1=re.split("\/"+tag,wd)[0]+'/'
queryrec=pd.read_csv(tmp1+"efetch_inputs/QueryToPMID_record.tsv",delimiter='\t')
queryrec['PMID'] = queryrec['PMID'].apply(str)
hypes['PMID'] = hypes['PMID'].apply(str)
hypes=hypes.join(queryrec.set_index('PMID'), on='PMID')
###
node_tbl=pd.read_csv(wd+"data/"+tag+"_allxall_nodes_table_bh_to_cytoscape.tsv",delimiter='\t') #"_node_table_with_rankings_strength_only.tsv",delimiter='\t')
#
hypes.to_csv(tmp+tag+"_Complete_edges_literature_links.tsv",sep='\t',index=False)
#
#
if 'incomplete' in sys.argv:
	#hypes.to_csv(tmp+tag+"_disc_Complete_edges_literature_links.tsv",sep='\t',index=False)
	node_tbl.to_csv(tmp+tag+"_disc_Complete_nodes_table_subgraph.tsv",sep='\t',index=False)
	edge_tbl[['kw1','kw2','pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration']].to_csv(tmp+tag+"_disc_Complete_edges_table_subgraph.tsv",sep='\t',index=False)
else:
	#hypes.to_csv(tmp+tag+"_Complete_edges_literature_links.tsv",sep='\t',index=False)
	node_tbl.to_csv(tmp+tag+"_Complete_nodes_table_subgraph.tsv",sep='\t',index=False)
	edge_tbl[['kw1','kw2','pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration']].to_csv(tmp+tag+"_Complete_edges_table_subgraph.tsv",sep='\t',index=False)

#node_tbl=pd.read_csv(wd+"data/"+tag+"_node_table_with_rankings.tsv",delimiter='\t')
cols=list(node_tbl.columns)
# category/label is 2nd column
cat=cols[1]

attr={node_tbl[cols[0]][i]:{cols[j]:node_tbl[cols[j]][i] for j in range(1,len(cols))} for i in range(len(node_tbl))}

# #converting node_tbl to arrtibutes dict
# attr={node_tbl['label'][i]:{
# 'category':node_tbl['category'][i],
# 'weight':node_tbl['node_weight'][i],
# 'cls_cnt':node_tbl['cls_cnt'][i],
# 'btw_cnt':node_tbl['btw_cnt'][i],
# 'strength':node_tbl['strength'][i]}
# #'Krank':node_tbl['Krank'][i],
# #'ranking':node_tbl['ranking'][i]}
# for i in range(len(node_tbl))}

# creating Network
Network=nx.from_pandas_edgelist(edge_tbl,'kw1','kw2',['pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration'])
#
gs=list(nx.connected_components(Network))
if len(gs) > 1:
	#print("Warning: Not every significant pair of entities is retained in the complete Network.\n(Probably some gene-only small subgraphs).")
	print("Warning: the network is not complete and so can possibly be genes and MeSH-only subgraphs")
	#gs={len(g):g for g in gs}
	#Network=Network.subgraph(gs[max(gs.keys())])
# adding node attributes
nx.set_node_attributes(Network,attr)
del attr
#
print('Network loaded')
#
def filter_nodes(x,fie,val):
    # x is
    node,att = x
    if att[fie] == val:
          return True  
    return False
### filter by label
meshs=list(filter(lambda node: filter_nodes(node,fie='category',val='MESH'),Network.nodes(data=True)))
genes=list(filter(lambda node: filter_nodes(node,fie='category',val='GENE'),Network.nodes(data=True)))
#
meshs=[t[0] for t in meshs]
genes=[t[0] for t in genes]
#
sub_meshs=Network.subgraph(meshs)
sub_genes=Network.subgraph(genes)
#
print('calculating centralities (Genes and Meshs networks)')
#
def netstats(g):
	closeness=nx.closeness_centrality(g,distance='inv_edge_weight')
	nx.set_node_attributes(g,closeness,"closeness_cnt")
	#
	betweenness=nx.betweenness_centrality(g,weight='inv_edge_weight')
	nx.set_node_attributes(g,betweenness,"betweenness_cnt")
	#
	strength=dict(g.degree(weight='R1C1'))
	nx.set_node_attributes(g,strength,'strength')
#
#netstats(sub_meshs)
#netstats(sub_genes)
#
### export gene-only and mesh-only
mesh_edge_tbl=nx.to_pandas_edgelist(sub_meshs,source='kw1',target='kw2')
mesh_node_tbl=pd.DataFrame.from_dict(dict(sub_meshs.nodes(data=True)), orient='index')
mesh_node_tbl=mesh_node_tbl.reset_index()
mesh_node_tbl=mesh_node_tbl.rename(columns={"index": "KW"})

mesh_node_tbl.to_csv(tmp+tag+"_MeSH_nodes_table_subgraph.tsv",sep='\t',index=False)
mesh_edge_tbl[['kw1','kw2','pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration']].to_csv(tmp+tag+"_MeSH_edges_table_subgraph.tsv",sep='\t',index=False)
#
gene_edge_tbl=nx.to_pandas_edgelist(sub_genes,source='kw1',target='kw2')
gene_node_tbl=pd.DataFrame.from_dict(dict(sub_genes.nodes(data=True)), orient='index')
gene_node_tbl=gene_node_tbl.reset_index()
gene_node_tbl=gene_node_tbl.rename(columns={"index": "KW"})
if 'incomplete' in sys.argv:
	gene_node_tbl.to_csv(tmp+tag+"_disc_Genes_nodes_table_subgraph.tsv",sep='\t',index=False)
	gene_edge_tbl[['kw1','kw2','pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration']].to_csv(tmp+tag+"_disc_Genes_edges_table_subgraph.tsv",sep='\t',index=False)

else:
	gene_node_tbl.to_csv(tmp+tag+"_Genes_nodes_table_subgraph.tsv",sep='\t',index=False)
	gene_edge_tbl[['kw1','kw2','pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration']].to_csv(tmp+tag+"_Genes_edges_table_subgraph.tsv",sep='\t',index=False)
#
edgeg=set(gene_edge_tbl.kw1.tolist()+gene_edge_tbl.kw2.tolist())
nodeg=set(gene_node_tbl.KW.tolist())
#
if len(edgeg) != len(nodeg):
	print("Warning: %i genes are disconnected from the gene-only network" % (len(nodeg - edgeg)))
	print(nodeg-edgeg)
	print("inspect possible connecting pathways...")
	inspect=list(nodeg-edgeg)
	inspect=['\n'+s for s in inspect]
	with open(wd+tag+"_Genes_unconnected.txt",'w') as file:
		file.writelines(inspect) 
#sys.exit("Gene-MeSH Network Completed - start connecting gene-only subgraphs")