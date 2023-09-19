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
from os import listdir
from os.path import isfile, join
from collections import Counter
from operator import itemgetter
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
thr= float(sys.argv[3])
### import Network after ranking nodes in R ###

edge_tbl=pd.read_csv(wd+"data/"+tag+"_allxall_sig_combinations_bh_to_cytoscape.tsv",delimiter='\t')
#ew_list=list(edge_tbl['edge_weight'])
edge_tbl['1-edge_weight']=[1-ew for ew in edge_tbl['edge_weight']]

node_tbl=pd.read_csv(wd+"data/"+tag+"_allxall_nodes_table_bh_to_cytoscape.tsv",delimiter='\t') #"_node_table_with_rankings_strength_only.tsv",delimiter='\t')

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
Network=nx.from_pandas_edgelist(edge_tbl,'kw1','kw2',['edge_weight','inv_edge_weight'])

# adding node attributes
nx.set_node_attributes(Network,attr)
del attr
#
print('Network loaded')
#
# ### finding motif GENE-GENE-MESH closed triplet ###
# def get_genelist(node,g,cat=cat):
# 	l=[]
# 	if node in g:
# 		for nb in g.neighbors(node):
# 			if g.nodes(cat)[nb]=='GENE':
# 				l.append(nb)
# 		if len(l)>=2:
# 			gene_nb=(node,l)
# 			return gene_nb

# # save as tuples gene_nb = (origin mesh,[list of genes])


# #
# def get_nodes(n1,n2,g=Network):
# 	if g.has_edge(n1,n2):
# 		return (n1,n2)

# # to return nodes instead of bool with starmap

# #
# def find_motif(mesh_gene_list_tuple,all_motifs):
# 	motifs=[]
# 	genes=mesh_gene_list_tuple[1]
# 	edge_list=[x for x in starmap(get_nodes,itertools.combinations(genes,2)) if x]
# 	if edge_list:
# 		for e in edge_list:
# 			motifs.append((e[0],e[1],mesh_gene_list_tuple[0]))
# 		if motifs:
# 			all_motifs.append(motifs)

# # returns motif as list of tuples(gene1,gene2,origin mesh)
# print("retrieve MeSH list")

# mesh_list=[node for node in list(Network.nodes) if Network.nodes(cat)[node]=='MESH']

# all_motifs=[]

# print("Find motifs")

# for mesh in mesh_list:
# 	g_genes=(get_genelist(mesh,Network))
# 	if g_genes:
# 		find_motif(g_genes,all_motifs)


# all_motifs=list(itertools.chain.from_iterable(all_motifs))

# print("Motifs extracted")

#all_motifs=pd.DataFrame(all_motifs, columns=['gene1','gene2','mesh'])

# LOAD SUBGRAPHS NODE TABLES TO ASSIGN SCORES AND SUBGRAPH BELONGING

# subdir=wd+'data/subgraphs/'#+tag+'/'
# nodetbls = [f for f in listdir(subdir) if bool(re.search("nodes", f))]
# nodetbls.sort()

# sub_nodes={}

# for f in nodetbls:
# 	ext=f.split('_')[-1]
# 	inn=int(ext.split('.')[0])
# 	f=pd.read_csv(subdir+f,delimiter='\t')
# 	sub_nodes[inn]=f

# ranked_mots=[]
# ser_kwpub=decompress_pickle(wd+"data/"+tag+"_KWtoPub_series.pbz2") # IT HAS TO BE UPDATED ALONG THE ITERATION
# tot_ids=len(set(itertools.chain.from_iterable(ser_kwpub.tolist())))
# #
# def represent(kw,tot_ids=tot_ids,ser=ser_kwpub):
# 	return len(ser[kw])/tot_ids
# #
# print("compute representation parameters...")
# #	
# def starrep(mot,sub,sub_nodes=sub_nodes):
# 	if all(m in sub_nodes[sub].iloc[:,0].tolist() for m in mot):
# 		# found the right subgraph
# 		subg=sub
# 		# now for the rankings
# 		df=sub_nodes[sub]
# 		#ranks=list(df.loc[df.iloc[:,0].isin(mot)].ranking)
# 		#score=hmean(ranks)
# 		# evaluate how "representative" the subgraph is of the "population" (articles) => PMIDs % of the total
# 		repGT=100*(len(set(itertools.chain.from_iterable(ser_kwpub[df.entity.tolist()].tolist())))/tot_ids)
# 		# evaluate, on average, how the nodes belonging to subgraph are "representative" of a "population" (articles) are => PMIDs % of the total
# 		repGM=100*(hmean(list(map(represent,df.entity.tolist()))))
# 		# evaluate how "representative" a motif is of the "population" (articles) => PMIDs % of the total
# 		repM=100*(len(set(itertools.chain.from_iterable(ser_kwpub[list(mot)].tolist())))/tot_ids)
# 		# harmonic mean motif-only
# 		repMM=100*(hmean(list(map(represent,ser_kwpub[list(mot)].index)))) 
# 		#
# 		return mot+(subg,repGT,repGM,repM,repMM)
# #
# ranked_mots=[i for i in starmap(starrep,itertools.product(all_motifs,list(sub_nodes.keys()))) if i]
# #
# # for mot,sub in itertools.product(all_motifs,list(sub_nodes.keys())):
# # 	if all(m in sub_nodes[sub].iloc[:,0].tolist() for m in mot):
# # 		# found the right subgraph
# # 		subg=sub
# # 		# now for the rankings
# # 		df=sub_nodes[sub]
# # 		ranks=list(df.loc[df.iloc[:,0].isin(mot)].ranking)
# # 		score=hmean(ranks)
# # 		# evaluate how "representative" the subgraph is of the "population" (articles) => PMIDs % of the total
# # 		repGT=100*(len(set(itertools.chain.from_iterable(ser_kwpub[df.label.tolist()].tolist())))/tot_ids)
# # 		# evaluate, on average, how the nodes belonging to subgraph are "representative" of a "population" (articles) are => PMIDs % of the total
# # 		repGM=100*(hmean(list(map(represent,df.label.tolist()))))
# # 		# evaluate how "representative" a motif is of the "population" (articles) => PMIDs % of the total
# # 		repM=100*(len(set(itertools.chain.from_iterable(ser_kwpub[list(mot)].tolist())))/tot_ids)
# # 		# harmonic mean motif-only
# # 		repMM=100*(hmean(list(map(represent,ser_kwpub[list(mot)].index)))) 

# # 		ranked_mots.append(mot+(subg,score,repGT,repGM,repM,repMM))

# #all_motifs=pd.DataFrame(ranked_mots, columns=['gene1','gene2','mesh','subgraph','ranking','rep_graphtot','rep_graphmean','rep_motif','rep_motif_mean'])
# all_motifs=pd.DataFrame(ranked_mots, columns=['gene1','gene2','mesh','subgraph','rep_graphtot','rep_graphmean','rep_motif','rep_motif_mean'])
# #
# print("...done")
# # rank by K-node score 

# # def add_scores(i,df=all_motifs):
# # 	score=sum([Network.nodes('Krank')[node] for node in df.loc[i]])
# # 	return score

# # all_motifs['score']=[add_scores(i) for i in range(len(all_motifs))]

# # all_motifs.sort_values(by='score',ascending=False,inplace=True)
# # all_motifs.reset_index(drop=True,inplace=True)

# # rank by ranking

# # def harm_mean(i, df=all_motifs):
# # 	ranks=[Network.nodes('ranking')[node] for node in df.loc[i]]
# # 	score=hmean(ranks)
# # 	return score


# # all_motifs['score']=[harm_mean(i) for i in range(len(all_motifs))]

# # FIND A WAY TO GET A PRE-DEFINED PRIORITY OF INTERSECTION 

# # define priority based on HARMONIC, AVERAGE REPRESENTATIVENESS (rep_GM)
# #lpri=list(set(all_motifs.ranking.tolist()))
# #lpri.sort()

# #all_motifs['priority']=all_motifs['ranking'].apply(lambda x: lpri.index(x))
# # sorty by (hierachically): priority, ranking, subgraph size
# all_motifs.sort_values(by=['rep_graphmean','subgraph'],ascending=True,inplace=True)
# all_motifs.reset_index(drop=True,inplace=True)
# all_motifs.to_csv(wd+"data/"+tag+"_motifs_subs_and_ranks.tsv",sep='\t')
# print('motifs saved')
# #
# ### STOP IF only 1 subgraph survives ###
# ### MOTIF REPRESENTATIVE THR ### 
# all_motifs=all_motifs[all_motifs['rep_graphtot']>=thr]
# all_motifs.reset_index(drop=True,inplace=True)
# #
# if len(set(all_motifs.subgraph.tolist())) < 2:
# 	print("There is just one subgraph having motifs and representing more than %i percent of papers" % (thr))
sd=re.split("tmp",wd)[0]
tmp=re.split(tag+"\/"+"data",wd)[0]
exec(open(sd+"code/sort_network.py").read())
sys.exit("Gene-MeSH Network Completed - start connecting gene-only subgraphs")
#print("Gene-MeSH Network Completed - start connecting gene-only subgraphs")
	#sys.exit("Gene-MeSH Network Completed - start connecting gene-only subgraphs")
###	
#all_motifs['subgraph']=random.choices([1,2,3],k=len(all_motifs))
###
# subgs=all_motifs.subgraph.tolist()
# subgs=list(set(subgs))
# subgs.sort()
# #
# indx_mot={s:list(all_motifs[all_motifs['subgraph']==s].index) for s in subgs}
# #
# query_ind=[]
# #
# # note: "second of the first is better than first of the second"
# # how many attemps for any pair of subgraphs?
# A=int(sys.argv[4])
# #
# def gapfind(s1,s2,attemps=A,im=indx_mot,mf=all_motifs):
# 	# create subsets
# 	#
# 	#
# 	s1df=mf.iloc[im[s1]] # list of indices
# 	s2df=mf.iloc[im[s2]]
# 	for s in [s1df,s2df]:
# 		s=s.sort_values('priority')
# 		s.reset_index(drop=True,inplace=True)
# 	# create priority queue 
# 	def pscore(i1,i2,s1=s1df,s2=s2df):
# 		#
# 		ind1=i1+1
# 		ind2=i2+1
# 		if min([ind1,ind2])==ind1:
# 			# it has priority 
# 			prio=0
# 		else:
# 			# (1,2) > (2,1)
# 			prio=1
# 		#
# 		k=hmean([ind1,ind2])
# 		return ((tuple(s1.loc[i1][:3]),tuple(s2.loc[i2][:3])),k,prio)
# 		#
# 	out=list(starmap(pscore,itertools.product(list(s1df.index)[:A],list(s2df.index)[:A])))
# 	out=sorted(out, key=itemgetter(1, 2))[:A]
# 	npairs=list(iter([t[0] for  t in out]))
# 	return {(s1,s2):npairs}	
# #
# query_ind=list(starmap(gapfind,itertools.combinations(subgs,2)))
# # #
# # i=0
# # while i < len(subgs):
# # 	s1=indx_mot[subgs[i]]
# # 	s2=dict(itertools.islice(indx_mot.items(), i+1,max(subgs)+1)).values() #indx_mot[subgs[i+1:]]
# # 	s2=list(itertools.chain.from_iterable(s2))
# # 	for j in s2:
# # 		query_ind+=list(itertools.product(s1,[j]))
# # 	i+=1

# # print("queries queue assessed")
# # def query_list(ind1,ind2,df=all_motifs):
# # 	m1=tuple(df.loc[ind1][:3])
# # 	m2=tuple(df.loc[ind2][:3])
# # 	return (m1,m2)

# # out=starmap(query_list,query_ind)
# compress_pickle(wd+'data/'+tag+'_query_queue',query_ind)