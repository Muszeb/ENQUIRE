#### WRAPPER SCRIPT TO EXECUTE ENQUIRE.sh ####
#### AND GENERATE AN OUTPUT REPRESENTATION IN THE FOR OF A PYTHON OBJECT ####

#### THESE VARIABLES SHALL BE PASSED TO ENQUIRE.sh (examples) ####
# wd=$(pwd)/
# tag=SplicingFactors_Neoplasms_Antigens
# ncores=32
# K=3
# A=2
# thr=2
# comb=4
# to_py=pmid-rnasplicin-antigens-set.txt
# etype='all'
# rscript=/home/musellla/miniconda3/envs/R363/bin/Rscript
####

import os
import sys
import numpy as np 
import pandas as pd
import subprocess
import glob 
import re
#
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
# Literature Data and all other edge/node tables #
def globtab(globbing):	
	litlinksfile=glob.glob(globbing)
	litlinksdata=[pd.read_csv(f,sep='\t') for f in litlinksfile]
	# regex names # 
	expind=[re.findall("^([^/]+)\\/(.*subgraphs_expansion)?([0-9]?[0-9]?)",s)[0][-1] for s in litlinksfile]
	expind=[int(x) if bool(re.search('[0-9]+',x)) else 0 for x in expind]
	return {k:v for k,v in zip(expind,litlinksdata)} 
#	
# subgraphs 
def globtabsub(globbing):	
	litlinksfile=glob.glob(globbing)
	litlinksdata=[pd.read_csv(f,sep='\t') for f in litlinksfile]
	# regex names # 
	expind=[re.findall("^([^/]+)\\/(.*subgraphs_expansion)?([0-9]?[0-9]?)",s)[0][-1] for s in litlinksfile]
	expind=[int(x) if bool(re.search('[0-9]+',x)) else 0 for x in expind]
	#
	subid=[int(re.findall("^.+_([0-9]+)\\.tsv$",s)[0][-1]) for s in litlinksfile]
	#
	tabl=pd.DataFrame({'Expansion':expind,'SubgraphId':subid,'File':litlinksfile})
	# nested table of subgraphs tables 
	return tabl.groupby('Expansion').apply(lambda x: {k:pd.read_csv(v,sep='\t') for k,v in zip(x.SubgraphId,x.File)}).to_dict() 


# 1) Pass Parameters to Python Variables (example) #

# wd='$(pwd)/'
# tag='SplicingFactors_Neoplasms_Antigens'
# ncores=32
# K=3
# A=2
# thr=2
# comb=4
# to_py='pmid-rnasplicin-antigens-set.txt'
# etype='all'
# rscript='/home/musellla/miniconda3/envs/R363/bin/Rscript'

# 2) Execute ENQUIRE.sh 

def run(tag,to_py,thr=1,wd="$(pwd)/",comb=4,A=2,K=3,etype='all',rscript="/home/musellla/miniconda3/envs/R363/bin/Rscript",ncores=32):
	#
	execstr="./code/ENQUIRE.sh -p %s -t %s -i %s -r %s -c %s -a %s -k %s -e %s -w %s -j %s"
	#
	process=subprocess.Popen('echo '+wd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
	wd,dnull=process.communicate()
	wd=wd.decode().strip('\n')
	#
	process=subprocess.Popen(execstr % (str(wd),
		str(tag),
		str(to_py),
		str(thr),
		str(comb),
		str(A),
		str(K),
		str(etype),
		str(rscript),
		str(ncores),), shell=True,
		stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	#
	stdout, stderr = process.communicate()
	# 3) Point at working directory (wd) and store dataframes into a dictionary #
	#
	if os.getcwd() != wd:
		os.chdir(wd+'tmp-'+tag)
	#
	ENQRES=dict()
	#
	# 4) save in dict-like structure #
	#
	ENQRES['Parameters']={'Working directory':wd,
		'Job description':tag,
		'Input corpus':to_py,
		'Representativeness (t)':thr,
		'Terms per query (k)':comb,
		'Maximum attempts (A)':A,
		'Expansion constraint (K)':K,
		'Entity types': etype,
		'Rscript':rscript,
		'Number of cores':ncores}
	#
	#
	ENQRES['Literature Data']=globtab('*/*Complete_edges_literature_links.tsv')
	ENQRES['Edgelist Data']=globtab('*/*Complete_edges_table_subgraph.tsv')
	ENQRES['Nodelist Data']=globtab('*/*Complete_nodes_table_subgraph.tsv')
	ENQRES['Subgraphs Node Data']=globtabsub('*/data/subgraphs/*nodes_table*')
	ENQRES['Subgraphs Edge Data']=globtabsub('*/data/subgraphs/*edges_table*')
	ENQRES['Gene-Subgraphs Edge Data']=globtabsub('*/data/gene-subgraphs/*edges_table*')
	ENQRES['Gene-Subgraphs Node Data']=globtabsub('*/data/gene-subgraphs/*nodes_table*')
	#
	return ENQRES
