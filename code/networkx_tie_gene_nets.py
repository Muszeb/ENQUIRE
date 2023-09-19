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
wd=sys.argv[1]
#
tag=sys.argv[2] #iter1
#
sd=re.split("tmp",wd)[0]
#tmp=re.split(tag+"\/",wd)[0]
tmp=re.split("\/"+tag,wd)[0]+'/'
#
prev=sys.argv[3] #"vcells"
prevwd=tmp+prev+"/"
#
sd=re.split("tmp",wd)[0]
tmp=''.join(re.split("\/"+tag,wd))
#
### load gene subgraphs ###
# dictionary of sub-gene-graphs (1 for edges 1 for nodes)
print("load sub-gene-graphs networks")
#
nod_fil=glob.glob(prevwd+'data/gene-subgraphs/*nodes*')
edg_fil=glob.glob(prevwd+'data/gene-subgraphs/*edges*')
#
nod_fil.sort()
edg_fil.sort()
#
subg_noded={i:pd.read_csv(j,delimiter='\t') for i,j in zip(range(0,len(nod_fil)),nod_fil)}
subg_edged={i:pd.read_csv(j,delimiter='\t') for i,j in zip(range(0,len(edg_fil)),edg_fil)}
#
def recnet(nod,edg):
	# nod and edg are dict iterms of the form (key,df)
	key,noddf = nod
	key,edgdf = edg
	cols=list(noddf.columns)
	if cols[0]=='label':
		cols[0]='KW'
	# category/label is 2nd column
	cols[1]='label'
	#
	noddf=noddf.rename(columns={n:c for n,c in zip(noddf.columns,cols)})
	#
	attr={noddf[cols[0]][i]:{cols[j]:noddf[cols[j]][i] for j in range(1,len(cols))} for i in range(len(noddf))}
	#
	net=nx.from_pandas_edgelist(edgdf,'kw1','kw2',['pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration'])
	# adding node attributes
	nx.set_node_attributes(net,attr)
	del attr
	#
	return (key,net)

#
#
subgnets=dict(starmap(recnet,zip(subg_noded.items(),subg_edged.items())))
#
### construct all possible pairs of suitable genes that would connect two subgraphs ### 
def gene_pairs(s1,s2):
	s1=list(s1.nodes())
	s2=list(s2.nodes())
	return list(itertools.product(s1,s2))
###
gp=list(starmap(gene_pairs,itertools.combinations(list(subgnets.values()),2)))	
gp=list(itertools.chain.from_iterable(gp))
#
preg=set(itertools.chain.from_iterable(gp))
#
print("all possible gene node pairs computed")
#
###
#
# load network resulted from paper expansion
#
print("load network resulted from papers' set expansion")
#
exp_edge_tbl=pd.read_csv(wd+"data/"+tag+"_allxall_sig_combinations_bh_to_cytoscape.tsv",sep='\t')
exp_node_tbl=pd.read_csv(wd+"data/"+tag+"_allxall_pre_nodes_table_bh_to_cytoscape.tsv",sep='\t')
exp_node_tbl=exp_node_tbl.rename(columns={'WEIGHT':'weight'})
#
cols=list(exp_node_tbl.columns)
# category/label is 2nd column
cat='label'
#
attr={exp_node_tbl[cols[0]][i]:{cols[j]:exp_node_tbl[cols[j]][i] for j in range(1,len(cols))} for i in range(len(exp_node_tbl))}
#
Network=nx.from_pandas_edgelist(exp_edge_tbl,'kw1','kw2',['OR','pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration'])
# adding node attributes
nx.set_node_attributes(Network,attr)
#del attr
#
genes=exp_node_tbl[exp_node_tbl['label']=='GENE']['KW'].tolist()
#
Net_genes=Network.subgraph(genes)
#
print("Gene Network and attributes loaded")
#
def spp(g1,g2,N=Network):
	if g1 in N.nodes() and g2 in N.nodes():
		try:
			return list(nx.all_shortest_paths(N,g1,g2))
		except:
			return None
#
print("computing connecting paths...")
shortsp=list(starmap(spp,gp))
shortsp=[p for p in shortsp if p]
shortsp=list(itertools.chain.from_iterable(shortsp))
addgenes=list(set(itertools.chain.from_iterable(shortsp)))

print("found %i new connecting genes" % (len(set(addgenes) - preg)))
#
# retrieve associated-only MeSH terms
#
meshs=exp_node_tbl[exp_node_tbl['label']=='MESH']['KW'].tolist()
#
def recmesh(g1,G=Network,M=meshs):
	neigh=list(nx.all_neighbors(G,g1))
	neigh=[m for m in neigh if m in M]
	return neigh
#
meshs=list(map(recmesh,addgenes))
meshs=list(set(itertools.chain.from_iterable(meshs)))
#
keepn=meshs+addgenes
#
keeNet=Network.subgraph(keepn)
#
### now create one big graph and update gene and mesh-only ones
print("combining networks...")
#
pre_edge_tbl=pd.read_csv(prevwd+prev+"_Complete_edges_table_subgraph.tsv",sep='\t')
pre_node_tbl=pd.read_csv(prevwd+prev+"_Complete_nodes_table_subgraph.tsv",sep='\t')
#
# the "exp_node_tbl" is behind in descriptiveness => adapt consequently pre node tbl
pre_node_tbl=pre_node_tbl[exp_node_tbl.columns]
#
cols=list(pre_node_tbl.columns)
# category/label is 2nd column
cat='label'
#
attr={pre_node_tbl[cols[0]][i]:{cols[j]:pre_node_tbl[cols[j]][i] for j in range(1,len(cols))} for i in range(len(pre_node_tbl))}
#
Network=nx.from_pandas_edgelist(pre_edge_tbl,'kw1','kw2',['OR','pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration'])
#
gs=list(nx.connected_components(Network))
if len(gs) > 1:
	print("Warning: Not every significant pair of entities is retained in the complete Network.\n(Probably some gene-only small subgraphs).")
	gs={len(g):g for g in gs}
	#unfrozen_graph = nx.Graph(frozen_graph)
	Network=nx.Graph(Network.subgraph(gs[max(gs.keys())]))
	#Network=Net
# adding node attributes
nx.set_node_attributes(Network,attr)
del attr
# new union graph
Network.update(edges=keeNet.edges(data=True),nodes=keeNet.nodes(data=True))
#
def netstats(g):
	closeness=nx.closeness_centrality(g,distance='inv_edge_weight')
	nx.set_node_attributes(g,closeness,"closeness_cnt")
	#
	betweenness=nx.betweenness_centrality(g,weight='inv_edge_weight')
	nx.set_node_attributes(g,betweenness,"betweenness_cnt")
	#
	strength=dict(g.degree(weight='edge_weight'))
	nx.set_node_attributes(g,strength,'strength')
#
print('calculating centralities (Complete Network)')
netstats(Network)
#
edge_tbl=nx.to_pandas_edgelist(Network,source='kw1',target='kw2')
node_tbl=pd.DataFrame.from_dict(dict(Network.nodes(data=True)), orient='index')
node_tbl=node_tbl.reset_index()
node_tbl=node_tbl.rename(columns={"index": "KW"})
#
### it MIGHT still be necessary to update merge_dict (ruling out unused papers) and
### recompute cross-weight, node-weight
dict_pub=decompress_pickle(tmp+"merge_dict_pub.pbz2")
#dict_not_pub=decompress_pickle(tmp+"merge_dict_not_pub.pbz2")
#
pmids=len(set(itertools.chain.from_iterable(list(dict_pub.values()))))
#
kws=node_tbl.KW.tolist()
pmin=[]
for k in kws:
	pmin+=dict_pub[k]
#
pmin=set(pmin) # sets of actually used papers
#
if len(pmin) < pmids:
	# update the dictionaries
	for k,val in dict_pub.items():
		dict_pub[k]=list(set(val) & pmin)


	#for k,val in dict_not_pub.items():
	#	dict_not_pub[k]=list(set(val) & pmin)
	#
	tups=[tuple(x) for x in edge_tbl[['kw1','kw2']].to_numpy()]
	#
	def inters(kw1,kw2,l,pmid=len(pmin),keys_pub=dict_pub):#,keys_not_pub=dict_not_pub):
		#	 TABLE APPEARENCE: #
		#		mesh2 	not-mesh2 #
		#						  #
		# mesh1 				  #
		# 						  #
		# not-mesh1				  #
		totalids=len(set(keys_pub[kw1]+keys_pub[kw2]))
		crossw=len((set(keys_pub[kw1]) & set(keys_pub[kw2])))/totalids # edge weight
		node1w=len(set(keys_pub[kw1]))/pmid
		node2w=len(set(keys_pub[kw2]))/pmid
		# save node weights info in a list of tuples
		l.append((kw1,node1w))
		l.append((kw2,node2w))
		#if mode=='edge_weight':
		return crossw		
			# return edge weight
		#else
	### 
	nodestab=[]
	edge_tbl['edge_weight'] = edge_tbl.apply(lambda x: inters(x['kw1'], x['kw2'],nodestab), axis=1)
	edge_tbl['inv_edge_weight']=[1-el for el in edge_tbl['edge_weight']]
	#
	#nodestab=list(set(itertools.chain.from_iterable(nodestab)))
	nodestab=pd.DataFrame(nodestab, columns=['KW','weight'])
	nodestab=nodestab.drop_duplicates()
#
# # additional mouse2human conversions and duplicate removal
# m2h_dict=decompress_pickle(sd+'input/mouse2human_genes_dict.pbz2')
# #
# ### keep track of m2h mappings (needed for co-sentence.py )
# m2h_map=[]
# ###
# def m2h_net(g1,m2h=m2h_dict,df=node_tbl.KW.tolist(),m2h_map=m2h_map):
# 	if g1 in m2h.keys():
# 		m2h_map.append((g1,m2h[g1]))
# 		return m2h[g1]
# 	elif g1.upper() in df:
# 		m2h_map.append((g1,g1.upper()))
# 		return g1.upper()
# 	else:
# 		m2h_map.append((g1,g1))
# 		return g1
# #
# node_tbl['KW']=node_tbl['KW'].apply(m2h_net)
# edge_tbl['kw1']=edge_tbl['kw1'].apply(m2h_net)
# edge_tbl['kw2']=edge_tbl['kw2'].apply(m2h_net)
# #
# m2h_map=dict({t[0]:t[1] for t in m2h_map})
# m2h_map=pd.Series(m2h_map)
# compress_pickle(tmp+"m2h_map_tie_nets",m2h_map)
# drop duplicates
# you have to consider swapped occurrences too
node_tbl=node_tbl.drop_duplicates(subset=['KW'])
#
oks=pd.DataFrame()
oks['kw1']=edge_tbl[['kw1','kw2']].apply(lambda x: x.sort_values()[0],axis=1)
oks['kw2']=edge_tbl[['kw1','kw2']].apply(lambda x: x.sort_values()[1],axis=1)
oks=list(oks.drop_duplicates().index)
edge_tbl=edge_tbl.filter(oks,axis=0)
#
node_tbl=node_tbl.reset_index(drop=True)
edge_node_tbl=edge_tbl.reset_index(drop=True)
#
node_tbl.to_csv(wd+tag+"_Complete_nodes_table_subgraph.tsv",sep='\t',index=False)
edge_tbl.to_csv(wd+tag+"_Complete_edges_table_subgraph.tsv",sep='\t',index=False)
# 
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
	return [{'PMID': i, 'kw1': v[0], 'kw2': v[1], 'link':s} for i,s in zip(l,links)]
###	
hypes=edge_tbl[['kw1','kw2']].apply(hyperlink,axis=1).tolist()
hypes=list(itertools.chain.from_iterable(hypes))
hypes=pd.DataFrame(hypes)
#
hypes.to_csv(wd+tag+"_Complete_edges_literature_links.tsv",sep='\t',index=False)
# recreate network
cols=list(node_tbl.columns)
# category/label is 2nd column
cat='label'
#
attr={node_tbl[cols[0]][i]:{cols[j]:node_tbl[cols[j]][i] for j in range(1,len(cols))} for i in range(len(node_tbl))}
#
Network=nx.from_pandas_edgelist(edge_tbl,'kw1','kw2',['OR','pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration'])
# adding node attributes
nx.set_node_attributes(Network,attr)
del attr
# genes- and MeSH-only
def filter_nodes(x,fie,val):
	node,att = x
	if att[fie] == val:
		return True  
	return False
#
meshs=list(filter(lambda node: filter_nodes(node,fie='label',val='MESH'),Network.nodes(data=True)))
genes=list(filter(lambda node: filter_nodes(node,fie='label',val='GENE'),Network.nodes(data=True)))
#
meshs=[t[0] for t in meshs]
genes=[t[0] for t in genes]
#
sub_meshs=Network.subgraph(meshs)
sub_genes=Network.subgraph(genes)
#
print('calculating centralities (Genes and Meshs networks)')
#
netstats(sub_meshs)
netstats(sub_genes)
#
### export gene-only and mesh-only
mesh_edge_tbl=nx.to_pandas_edgelist(sub_meshs,source='kw1',target='kw2')
mesh_node_tbl=pd.DataFrame.from_dict(dict(sub_meshs.nodes(data=True)), orient='index')
mesh_node_tbl=mesh_node_tbl.reset_index()
mesh_node_tbl=mesh_node_tbl.rename(columns={"index": "KW"})
mesh_node_tbl.to_csv(wd+tag+"_MeSH_nodes_table_subgraph.tsv",sep='\t',index=False)
mesh_edge_tbl.to_csv(wd+tag+"_MeSH_edges_table_subgraph.tsv",sep='\t',index=False)
#
gene_edge_tbl=nx.to_pandas_edgelist(sub_genes,source='kw1',target='kw2')
gene_node_tbl=pd.DataFrame.from_dict(dict(sub_genes.nodes(data=True)), orient='index')
gene_node_tbl=gene_node_tbl.reset_index()
gene_node_tbl=gene_node_tbl.rename(columns={"index": "KW"})
gene_node_tbl.to_csv(wd+tag+"_Genes_nodes_table_subgraph.tsv",sep='\t',index=False)
gene_edge_tbl.to_csv(wd+tag+"_Genes_edges_table_subgraph.tsv",sep='\t',index=False)
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
#
# def lastconnect(g,N=Network,gs=edgeg,A=3):
# 	it=itertools.product([N],[g],list(gs))
# 	out=starmap(nx.all_shortest_paths,it)
# 	out=[list(o) for o in out]
# 	out=list(itertools.chain.from_iterable(out))
# 	out.sort(key=len)
# 	# select connections based on A
# 	i=0
# 	j=0
# 	taken=[]
# 	res=[]
# 	while i < A and j < len(out):
# 		if out[j][-1] not in taken:
# 			taken.append(out[j][-1])
# 			res.append(out[j])
# 			i+=1
# 			j+=1
# 		else:
# 			j+=1
# 	#
# 	return(res)

# out=map(lastconnect,nodeg-edgeg)
