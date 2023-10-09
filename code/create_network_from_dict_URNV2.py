# CREATE NETWORK FROM DICT
# exec(open("/home/musellla/tam_textmining/code/create_network_from_dict.py").read())
import os
import sys
import numpy as np 
import pandas as pd
import itertools
from itertools import islice,product,starmap
import multiprocessing as mp
from multiprocessing import Pool
from collections import OrderedDict
import re
import random
from random import sample
import time
import networkx as nx
import networkx.algorithms.community as nx_comm
import bz2
import json
from multiprocess import Process, Manager
from scipy.stats import fisher_exact
from scipy.stats import rankdata
from scipy.stats import norm,poisson
import scipy.special
import threading
from time import sleep
from datetime import datetime
import scipy
import csv
from datetime import datetime
from scipy.stats import binom
from scipy.stats import hypergeom
import copy
import statsmodels.api as sm
import scipy.special
import math
from itertools import compress
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
sd= re.split("tmp",wd)[0]
tag=sys.argv[2] #"iter1"
hd = re.split("\/"+tag,wd)[0]+'/'
ncores= int(sys.argv[3]) #32 #round(lim/tot)

#### FISHER TEST #####

# load dictionaries kw:list of pmis

keys_pub=decompress_pickle(wd+'data/'+tag+'_ids_to_kw.pbz2') 

print('mesh-terms loaded')

if len(keys_pub)==0:
	print("WARNING: MeSH dictionary is empty") # context pruning based on mesh tree
else:
	mesh_dict=decompress_pickle(sd+'input/meshD_names_to_ids.pbz2')

	# reference translation
	ref_mesh={key.replace(',','').lower():key for key in mesh_dict.keys()}
	#
	exclude=['B','E','K','F','G01','G02','G17','H','N','I','L','V','J','Z','M'] # uninteresting mesh categories


	#remove commas
	keylist=list(mesh_dict.keys())
	mesh_dict={key.replace(',','').lower():mesh_dict[key] for key in keylist}

	prots=decompress_pickle(sd+'input/prots-MESH_to_prune.pbz2')
	prots=[p.lower() for p in prots]


	# prune motives with uninteresting contexts (see exclude list)
	def prune_con(mot,mesh_dict,ex,prots):
		mesh=mot
		try:
			ids=mesh_dict[mesh.split('/')[0]] # it's a list
			# exclude predefined categories
			def start_ex(idd,exw):
				return idd.startswith(exw)
			# this way
			# we can therefore exclude any prefix
			combs=itertools.product(ids,ex)
			out=starmap(start_ex,combs)
			#
			if any(out):
				return None
			# exclude protein names
			elif mesh in prots:
				return None
			# otherwise pass
			else:
				return mot
		except:
			print("warning: not a MeSH term: ", mesh)
			return None
	# remove commas from keys_pub
	mesh_keys=list(keys_pub.keys())
	keys_pub={key.replace(',','').lower():keys_pub[key] for key in mesh_keys}

	mesh_keys=list(keys_pub.keys())
	new_keys=[prune_con(k,mesh_dict,exclude,prots) for k in mesh_keys]
	new_keys=[el for el in new_keys if el]

	d={}

	for k in range(len(new_keys)):
		if new_keys[k] in keys_pub:
			d[new_keys[k]]=keys_pub[new_keys[k]]


	keys_pub=d

	del d
	print('pruned mesh terms')

#load genes and convert mouse to human

genes_pub=decompress_pickle(wd+"data/"+tag+"_ids_to_genes.pbz2")
#
if len(genes_pub)==0:
	print("WARNING: gene dictionary is empty")
else:
	#	
	m2h_dict=decompress_pickle(sd+'input/mouse2human_genes_dict.pbz2')

	gene_keys=list(genes_pub.keys())

	print(len([g for g in gene_keys if g.isupper()]), "uppercase genes before m2h transform")

	for key in gene_keys:
		pm=genes_pub[key]
		if key.isupper():
			pass
		elif key in m2h_dict:
			#pm=genes_pub[key]
			if m2h_dict[key] in genes_pub: 
				hum=m2h_dict[key]
				genes_pub[hum]+=pm
				#genes_pub[hum]=list(set(genes_pub[hum]))
				del genes_pub[key]
		#or simply as uppercase:	
		elif key.upper() in genes_pub and key.upper() != key: # OTHERWISE YOU SIMPLY DELETE THE ENTRY!
			hum=key.upper()
			genes_pub[hum]+=pm
			#genes_pub[hum]=list(set(genes_pub[hum]))
			del genes_pub[key]

	print('genes loaded')

	print(len([g for g in list(genes_pub.keys()) if g.isupper()]), "uppercase genes after m2h transform")

	#del m2h
	del m2h_dict
	del gene_keys

#create list of all pmids
idvalues=list(keys_pub.values())+list(genes_pub.values())
pmids=list(itertools.chain.from_iterable(idvalues))
pmids=list(set(pmids))

# # create the "negation of it"
# keys_not_pub=dict({key:[] for key in keys_pub.keys()})
# for key in keys_not_pub:
# 	keys_not_pub[key]= list(set(pmids) - set(keys_pub[key]))

# genes_not_pub=dict({key:[] for key in genes_pub.keys()})
# for key in genes_not_pub:
# 	genes_not_pub[key]= list(set(pmids) - set(genes_pub[key]))
#

#merge kw and gene dictonary save lables to apply to nodes later

labels=[("MESH",k)for k in keys_pub.keys()]+[("GENE", g) for g in genes_pub.keys()] 


merge_dict_pub = {**keys_pub, **genes_pub}
#merge_dict_not_pub = {**keys_not_pub, **genes_not_pub}

# save occurrence dict as series
# INCLUDE previous iteration's occurrence dict (if it>1)
# save lables to apply to nodes laterit=sys.argv[4]
it=sys.argv[4]

if int(it) > 0: 
	
	print("merge iterations")
	# with open(hd+"previous_iteration.txt","r") as file:
	# 	prev=file.readline()
	# 	prev=prev.split('\n')[0]
	
	prev_dict_pub=decompress_pickle(hd+"merge_dict_pub.pbz2")
	#prev_dict_not_pub=decompress_pickle(hd+"merge_dict_not_pub.pbz2")
	keyl_new_pub=list(merge_dict_pub.keys())
	#keyl_new_not_pub=list(merge_dict_not_pub.keys())
	prev_labels=decompress_pickle(hd+"kw_labels.pbz2")
	labels+=prev_labels
	labels=list(set(labels))
	for key in keyl_new_pub:
		if key in prev_dict_pub:
			prev_dict_pub[key]=list(set(prev_dict_pub[key]+merge_dict_pub[key]))
		else:
			prev_dict_pub[key]=merge_dict_pub[key]
	
	# for key in keyl_new_not_pub:
	# 	if key in prev_dict_not_pub:
	# 		prev_dict_not_pub[key]=list(set(prev_dict_not_pub[key]+merge_dict_not_pub[key]))
	# 	else:
	# 		prev_dict_not_pub[key]=merge_dict_not_pub[key]		
	
	merge_dict_pub=prev_dict_pub
	#merge_dict_not_pub=prev_dict_not_pub
	
	# if under the "gene subgraphs expansion", neglects tests on MeSH terms ?
	#
	# # update list of all pmids
	print('dictionaries merged')		
#
#pmids=list(merge_dict_pub.keys())
#
compress_pickle(hd+"kw_labels", labels)
compress_pickle(hd+"merge_dict_pub",merge_dict_pub)
#compress_pickle(hd+"merge_dict_not_pub",merge_dict_not_pub)
compress_pickle(wd+"data/"+tag+"_KWtoPub_series",pd.Series(merge_dict_pub))
bu_merge_dict_pub = copy.deepcopy(merge_dict_pub)
#bu_merge_dict_not_pub = copy.deepcopy(merge_dict_not_pub)
print('dictionaries saved')
### delete the unmerged dicts

del keys_pub
#del keys_not_pub
del genes_pub
#del genes_not_pub
##### CREATE CO-OCC MATRIX #####
test_keys=list(merge_dict_pub.keys()) # test only newly found entities
genes=list(merge_dict_pub.keys())
pmids=list(set(itertools.chain.from_iterable(merge_dict_pub.values())))
#### CONFIGURATION MODEL AS AN URN PROBLEM ####
### MENTIONS CANNOT BE INTERPRETED AS NODE DEGREES! IT'S A VERTEX-CENTRIC MEASUREMENT
print("annotate co-occurrences...")
#
var_names_in=["dict_pub_"+str(i) for i in range(ncores)]

for name in var_names_in:
	globals()[name] = merge_dict_pub

combs=list(itertools.combinations(test_keys,2))

print("%i combinations created" % (len(combs),))

# create chunks of size ncores
print("create edge_list...")
#chunks = iter([test_keys[x:x+ncores] for x in range(0, len(test_keys), ncores)]) #iter([combs[x:x+ncores] for x in range(0, 10000, ncores)])

def edge_list(kw1,kw2,keys_pub,N=len(pmids),path=wd+"data/"+tag+"_edge_list_allxall.tsv"):
	#
	#kw1,kw2=kws
	keys_pub=globals()[keys_pub]	
	#
	global_lock = threading.Lock()	
	#
	R=len(keys_pub[kw1]) # N1
	C=len(keys_pub[kw2]) # N2
	#
	cross=set(keys_pub[kw1]) & set(keys_pub[kw2])
	#
	R1C1=len(cross) # J
	#print(R1C1)
	#print("R1C1:",R1C1)
	R1C2=R - R1C1 #len((set(keys_pub[kw1]) & set(keys_not_pub[kw2])))
	#print("R1C2:",R1C2)
	R2C1=C - R1C1 #len((set(keys_not_pub[kw1]) & set(keys_pub[kw2])))
	#print("R2C1",R2C1)		
	R2C2=N - C - R1C2 #len((set(keys_not_pub[kw1]) & set(keys_not_pub[kw2])))
	#print("R2C2:",R2C2)				
	table = np.array([[R1C1, R1C2], [R2C1, R2C2]])
	#
	# if any((x <= 0 for x in [R1C1, R1C2, R2C1, R2C2])):
	# 	print(N,R,C,R1C1,table)
	if any((x <= 0 for x in [(R1C1+R2C1),(R1C2+R2C2),(R1C1+R1C2),(R2C1+R2C2)])):
		edge_weight=0
	else:
		#P_ab=R1C1/(R1C1+R2C1) - R1C2/(R1C2+R2C2)
		#P_ba=R1C1/(R1C1+R1C2) - R2C1/(R2C1+R2C2)
		edge_weight=fisher_exact(table, alternative='greater').pvalue # P_ba*P_ab
	#
	#pg=np.prod([urn(N,E,Ka,Ci)*urn(N,E,Kb,Ci) for Ci in Cis])
	# #return (gene,kw,oddsr,p2s,pl,pg)
	#if oddsr != 0 and oddsr != float('inf'):
	out=iter([(kw1,kw2,edge_weight,pmid) for pmid in cross])
	#csvfile.writerow(res)
	#d[(kw1,kw2)]=res
	#	
	#out=iter([gen for gen in starmap(contin,combs) if gen])
	def writerows(self, rowdicts):
		return self.writer.writerows(map(self._dict_to_list, rowdicts))
	with global_lock:			
		with open(path,"a") as file:
			#print(datetime.now(),"I am writing edges for:",kw1,'AND',kw2)
			csvfile=csv.writer(file,delimiter='\t')
			csvfile.writerows(out)
#
def distro(datal,ncores,d1):
	# d1 and d2 are strings
	# set repetitions
	rep=int(len(datal)/ncores)+2
	d1s=[d1+str(i) for i in range(0,ncores)]*rep
	#d2s=[d2+str(i) for i in range(0,ncores)]*rep
	#d3s=[d3+str(i) for i in range(0,ncores)]*rep
	#print(d1s)
	#d2s=[d2+str(i) for i range(0,ncores)]*rep
	return iter((datal[i][0],datal[i][1],d1s[i]) for i in range(0,len(datal)))
#
start=time.time()
with Pool(ncores) as pool:
	#res=[pool.apply(find_in_series,args=(idd,globals()["abs_dict_"+str(ids.index(idd))],globals()["syn_dict_"+str(ids.index(idd))])) for idd in ids]   
	print("starmap...")
	abs_dict=pool.starmap(edge_list,distro(combs,ncores,'dict_pub_'))
	print("... edgelist stored in TSV file")
	#res=dict(res)
#
del combs
print("completed in %2f seconds" % (time.time()-start))
#
print("load multi graph...")
#edges=pd.read_csv(wd+"data/"+tag+"_edge_list_allxall.tsv",delimiter='\t',header=None,names=["kw1","kw2","edge_weight",'PMID'])
try:
	multig=nx.read_edgelist(wd+"data/"+tag+"_edge_list_allxall.tsv",delimiter='\t', create_using=nx.MultiGraph, nodetype=str, data=(("edge_weight", float),("PMID",int)), edgetype=None, encoding='utf-8')
except FileNotFoundError:
	print("STOP: there is no edgelist to construct the multigraph")
	print("...the multi graph has no edges, exiting...")
	sys.exit(11)
#abort if there are no edges
print("checking edges...")
if multig.size()==0:
	print("...the multi graph has no edges, exiting...")
	sys.exit(11)
else:
	print("...passed")
	print("graph size:", multig.size())

#
### 
print("compute degrees, csi matrix, M, m...")
#
### HENCE EVERY PMID MUST CONTAIN AT LEAST 2 ENTITIES BEING MENTIONED 
degseq=dict(multig.degree())
degs=[d for k,d in degseq.items()]  # degree sequence
vseq=pd.Series({v:i for v,i in zip(degseq.keys(),range(len(degseq)))})
### Number of edges is the sum of the degree sequence divided by two
m=int(sum(degs)/2) # or multig.size(), but in this way we check that everything is alright #
### COMPUTE CSI MATRIX, WHICH EXPLAINS ALL POSSIBLE COMBINATIONS GIVEN TWO NODES DEGREE ###
csi=np.outer(np.array(degs), np.array(degs))
### NO SELF LOOPS, REMOVE DIAGONAL VALUES 
np.fill_diagonal(csi,0)
### M is the sum of all realizations of stubs combinations (in any order)
M=np.sum(csi)
### Only thing left is the adjacency matrix - Aij can be accessed on the fly  
A=pd.Series(list(multig.edges())).value_counts().to_dict()
###
### HYPERGEOM STAT ###
def hyperurn(Aij,csi,vseq,M=M,m=m):
	csi=globals()[csi]
	vseq=globals()[vseq]
	# Aij is a dict item, annotated as a ((v1,v2),count)
	ij,Aij=Aij
	poix=vseq[list(ij)]
	#csiij=2*csi[poix[0],][poix[1]]
	csiij=2*csi[poix.iloc[0]][poix.iloc[1]]
	return ij,1-hypergeom.cdf(Aij, M, csiij, m)

###
### CONTINOUS APPROXIMATION ###
# Let X ∼ Hypergeometric ⁡ ( N , K , n ) and p=K/N. 
#     If n n is large, N N and K K are large compared to n n, and p p is not close to 0 or 1, then

#         P ( X ≤ k ) ≈ Φ ( k − n p n p ( 1 − p ) ) P(X\leq k)\approx \Phi \left({\frac {k-np}{\sqrt {np(1-p)}}}\right)

# where Φ \Phi is the standard normal distribution function 
#
def approxurn(Aij,csi,vseq,M=M,m=m):
	csi=globals()[csi]
	vseq=globals()[vseq]
	# Aij is a dict item, annotated as a ((v1,v2),count)
	ij,Aij=Aij
	poix=vseq[list(ij)]
	#csiij=2*csi[poix[0],][poix[1]]
	csiij=2*csi[poix.iloc[0]][poix.iloc[1]]
	p=csiij/M
	Z=(Aij-m*p)/math.sqrt(((M-m)/(M-1))*m*p*(1-p))
	return ij,1-norm.cdf(Z)

####
var_names_in=["csi_"+str(i) for i in range(ncores)]

for name in var_names_in:
	globals()[name] = csi

var_names_in=["vseq_"+str(i) for i in range(ncores)]

for name in var_names_in:
	globals()[name] = vseq	

def distro(datal,ncores,d1,d2):
	# d1 and d2 are strings
	# set repetitions
	rep=int(len(datal)/ncores)+2
	d1s=[d1+str(i) for i in range(0,ncores)]*rep
	d2s=[d2+str(i) for i in range(0,ncores)]*rep
	#d3s=[d3+str(i) for i in range(0,ncores)]*rep
	#print(d1s)
	#d2s=[d2+str(i) for i range(0,ncores)]*rep
	return iter((datal[i],d1s[i],d2s[i]) for i in range(0,len(datal)))

print("compute urn/configuration model statistics")


if math.isinf(scipy.special.comb(M,m)):
	print("WARNING: cannot compute exact statistics\nUsing normal approximation...")	
	start=time.time()
	with Pool(ncores) as pool:
		#res=[pool.apply(find_in_series,args=(idd,globals()["abs_dict_"+str(ids.index(idd))],globals()["syn_dict_"+str(ids.index(idd))])) for idd in ids]   
		print("starmap...")
		stats=pool.starmap(approxurn,distro(list(A.items()),ncores,'csi_','vseq_'))
		print("... done in %2f seconds" % (time.time()-start))
	#res=dict(res)
else:
	start=time.time()
	with Pool(ncores) as pool:
		#res=[pool.apply(find_in_series,args=(idd,globals()["abs_dict_"+str(ids.index(idd))],globals()["syn_dict_"+str(ids.index(idd))])) for idd in ids]   
		print("starmap...")
		stats=pool.starmap(hyperurn,distro(list(A.items()),ncores,'csi_','vseq_'))
		print("... done in %2f seconds" % (time.time()-start))
#
# ADJUST PVALUES
print("pvalue adjustment")
#
qvals=sm.stats.multipletests([s[1] for s in stats],method='fdr_bh',alpha=0.01)[1].tolist()
qvals=[(t[0],q) for t,q in zip(stats,qvals)]
#qvals2=sm.stats.multipletests([s[1] for s in stats],method='holm-sidak')[1]
#qvals2=[(t[0],q) for t,q in zip(stats,qvals2)]
### Assign p- and q- values to the network ###
### for multigraphs, we need a (v1,v2,key) tuple as dict keys for the attributes to assign ###
def multiatt(stat,d):
	ij,stat=stat
	K=A[ij]
	i,j=ij
	for k in range(K):
		d[(i,j,k)]=stat 
#
pdict=dict()
x=[multiatt(stat,pdict) for stat in stats]
#
qdict=dict()
x=[multiatt(qval,qdict) for qval in qvals]
#
# add also R1C1 info 
adict=dict()
x=[multiatt(aij,adict) for aij in [(k,v) for k,v in A.items()]]
## assign!
nx.set_edge_attributes(multig, pdict, "pval")
nx.set_edge_attributes(multig, qdict, "qval")
nx.set_edge_attributes(multig, adict, "R1C1")
## compute Z score for R1C1 and use 1-Phi(Z) as edge weight ## 
mR1C1=np.mean(list(nx.get_edge_attributes(multig,'R1C1').values()))
#sdR1C1=np.std(list(nx.get_edge_attributes(multig,'R1C1').values()))
## inv_edge_weight as e^(-R1C1+1) ##
def multiatt(stat,d):
	ij,stat=stat
	K=A[ij]
	i,j=ij
	for k in range(K):
		#d[(i,j,k)]=math.exp(1-stat)
		# Zero Truncated Poisson: P( X>=k | X>0 )
		d[(i,j,k)]=(1-poisson.cdf(stat,mR1C1))/(1-poisson.cdf(0,mR1C1)) 
#
idict=dict()
x=[multiatt(aij,idict) for aij in [(k,v) for k,v in A.items()]]
nx.set_edge_attributes(multig, idict, "inv_edge_weight")
# used for louvain community detections
idict={k:1-v for k,v in idict.items()}
nx.set_edge_attributes(multig, idict, "weight")
## todo: compund PMIDs column
#dict(nx.get_edge_attributes(multig, "PMID"))
##
print('pruning and creating edge table...')
### PRUNE BY QVALUE
bade=[k for k,v in qdict.items() if v>=0.01]
#bade1=[k for k,v in adict.items() if v<3]
#bade=list(set(bade).union(set(bade1)))
prunedg=copy.deepcopy(multig)
print("... by adjusted p-values...")
prunedg.remove_edges_from(bade)
prunedg=nx.Graph(prunedg)
#
#sdR1C1=np.std(list(nx.get_edge_attributes(prunedg,'R1C1').values()))
#
# def pois(param,sample_data=np.array(list(nx.get_edge_attributes(prunedg,'R1C1').values()))):
#     # Calculate negative log likelihood
# 	nll = -np.sum(poisson.logpmf(sample_data, mu=param))#, scale=sd))
# 	return nll
# #
# results = minimize(pois, mR1C1, method='Nelder-Mead')
## inv_edge_weight as e^(-R1C1+1) ##
#
if  prunedg.size()!=0:
	print("...and by deleting Louvain-communities of size 1 ...")
	mR1C1=np.mean(list(nx.get_edge_attributes(multig,'R1C1').values()))
	def multiatt(stat,d):
		ij,stat=stat
		K=A[ij]
		i,j=ij
		#for k in range(K):
		#d[(i,j)]=math.exp(1-stat)
		#d[(i,j)]=1-norm.cdf(stat-mR1C1)/sdR1C1 
		# Zero Truncated Poisson: P( X>=k | X>0 )		
		d[(i,j)]=(1-poisson.cdf(stat,mR1C1))/(1-poisson.cdf(0,mR1C1)) 
	#
	idict=dict()
	x=[multiatt(aij,idict) for aij in [(k,v) for k,v in A.items()]]
	nx.set_edge_attributes(prunedg, idict, "inv_edge_weight")
	idict={k:-1*math.log10(max([v,1e-15])) for k,v in idict.items()}
	#idict={k:1-v for k,v in idict.items()}
	nx.set_edge_attributes(prunedg, idict, "weight")
	try:
		random.seed(2202)
		np.random.seed(2202)                                                         
      	os.environ['PYTHONHASHSEED'] = str(2202)
		CC=nx_comm.louvain_communities(prunedg,weight='weight',seed=2202)
	except ZeroDivisionError:
		random.seed(2202)
		np.random.seed(2202)                                                         
      	os.environ['PYTHONHASHSEED'] = str(2202)
		print("WARNING!\nAll edges have a weight of 0, hence a weighted Louvain clustering is not possible.")
		CC=nx_comm.louvain_communities(prunedg,weight=None,seed=2202)
	badv=[c for c in CC if len(c)==1]
	if len(badv)>0:
		badv=set.union(*[c for c in CC if len(c)==1])
		print("found %i communities and %i badly-connected nodes" % (len(CC),len(badv)))
		prunedg.remove_nodes_from(badv)
	else:
		print("found %i communities and %i badly-connected nodes" % (len(CC),len(badv)))
	print("...pruning done")
else:
	print("...after pruning there are no edges left, exiting...")
	print("graph size: ",prunedg.size())
	sys.exit(11)
#

res = nx.to_pandas_edgelist(prunedg)
### COUNTER BALANCE TEST STATISTICS: e^-(R1C1-1) ###
#res['inv_edge_weight']=res.R1C1.apply(lambda x: math.exp(-x+1))
### RENAME COLUMNS 
res.rename({'source': 'kw1','target': 'kw2','pval':'pg','qval':'crit_p'},axis = "columns", inplace = True)  

##################################################################

##### CREATE EDGES AND NODE ATTRIBUTES #####

tups=[tuple(x) for x in res[['kw1','kw2']].to_numpy()]
def inters(kw1,kw2,l,pmid=len(pmids),keys_pub=bu_merge_dict_pub):#keys_not_pub=bu_merge_dict_not_pub):
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
	l.append((kw1,node1w,len(set(keys_pub[kw1]))))
	l.append((kw2,node2w,len(set(keys_pub[kw2]))))
	#if mode=='edge_weight':
	return crossw		
		# return edge weight
	#else
### 
nodestab=[]
x = res.apply(lambda x: inters(x['kw1'], x['kw2'],nodestab), axis=1) # res['edge_weight']
#
#nodestab=list(set(itertools.chain.from_iterable(nodestab)))
nodestab=pd.DataFrame(nodestab, columns=['KW','node_weight','occurrence'])
nodestab=nodestab.drop_duplicates()

#add labels
labels=dict({t[1]:t[0] for t in labels})

kw_want=list(nodestab["KW"])

li_lab=[]

for kw in kw_want:
	li_lab.append(labels[kw])

nodestab["category"]=li_lab
nodestab=nodestab[['KW','category','node_weight','occurrence']] # reorder for consistency
#
# import edge and node tbl
edge_tbl=res #pd.read_csv(wd+"data/"+tag+"allxall_sig_combinations_bh_to_cytoscape.tsv",delimiter='\t')
node_tbl=nodestab #pd.read_csv(wd+"data/"+tag+"allxall_nodes_table_bh_to_cytoscape.tsv",delimiter='\t')
print(edge_tbl.head())
#edge_tbl['inv_edge_weight']= [-math.log10(el) for el in edge_tbl['edge_weight']] #[1-el for el in edge_tbl['edge_weight']]
edge_tbl=edge_tbl.drop_duplicates(subset=['kw1','kw2'])
# add iteration information
edge_tbl['iteration']=[int(it)] * len(edge_tbl)
#
edge_tbl=edge_tbl[["kw1","kw2",'pg','crit_p','edge_weight','inv_edge_weight','R1C1','iteration']]
#
print("save edge table")
edge_tbl.to_csv(wd+"data/"+tag+"_allxall_sig_combinations_bh_to_cytoscape.tsv",sep='\t',index=False)
node_tbl.to_csv(wd+"data/"+tag+"_allxall_pre_nodes_table_bh_to_cytoscape.tsv",sep='\t',index=False)
###