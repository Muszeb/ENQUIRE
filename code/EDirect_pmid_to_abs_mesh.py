#!/home/musellla/miniconda3/envs/wokenv/bin/python
# PMIDS to Abstract MESH
# exec(open("/home/musellla/tam_textmining/code/pmid_to_abs_mesh.py").read())
import os
import sys
#from nltk.corpus.reader.plaintext import PlaintextCorpusReader
import numpy as np 
import pandas as pd
import itertools
import multiprocessing as mp
from collections import OrderedDict
from nltk.stem import WordNetLemmatizer
import re
import random
from Bio import Entrez
from Bio.Entrez import efetch,esearch, read
from scipy.spatial import distance_matrix
import networkx as nx
import bz2
import json
import gzip
from multiprocess import Pool
import time
from itertools import islice
import nltk
nltk.download('wordnet')
nltk.download('stopwords')
from nltk.corpus import stopwords
import string
nltk.download('words')
from nltk.corpus import words
from nltk.stem import WordNetLemmatizer
from english_words import english_words_set
import spacy
import scispacy
import xmltodict
import dpath
from collections import ChainMap
import random
#import abbreviations
#from abbreviations import schwartz_hearst
#import nltk
#from sklearn.datasets import load_files
#nltk.download('wordnet')
#nltk.download('stopwords')
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


wd= sys.argv[1] #'/home/musellla/tam_textmining/'
tag= sys.argv[2]
idlist= sys.argv[3]
ncores= int(sys.argv[4])
#
sd=re.split("tmp",wd)[0]
#
### HANDLE ENTITY-TYPE RESTRICTION OPTION ### 
try:
    etype=sys.argv[5]
except Exception: 
    pass

#load batch of ids (a .txt file)
# decompress_pickle(wd+'data/pmid_queries_to_efetch.pbz2')
#idlist=idlist[:20]
with open(idlist,"r") as file:
	idlist=file.read().splitlines() #file.readlines()	
#
try:
	idlist=[int(i) for i in idlist]
except:
	sys.exit("ERROR: INVALID INPUT FILE")
#
print("list of PMIDs loaded")
#
id_abs={}
#
#ip(tam_id,tam_tx)}
print("fetch PMIDs")
#
def fetchlist(id2q,m):
	#keysform=["'"+v+"'" for v in keysform]
	id2q=[str(i) for i in id2q]
	print("try posting and fetching %i PMIDs..." % (len(id2q))) # .split(','))))
	if len(id2q)<m:		
		#
		id2q=','.join([str(i) for i in id2q])
		# efetch 1st kw, efilter with the others, then filter out reviews #
		query="epost -db pubmed -id '%s' | efetch -format xml > %s" % (id2q, wd+'data/'+tag+'_Efetch_Results.xml')#" efilter -query " % (webenvfile,keysform[0])
		os.system(query)
		with open(wd+'data/'+tag+'_Efetch_Results.xml') as xml:
			return dict({0:xmltodict.parse(xml.read())})
	else:
		print("generate %i PMIDs-long chunks..." % (m))
		chunks=[','.join(id2q[x:x+m]) for x in range(0, len(id2q), m)]
		queries=["epost -db pubmed -id '%s' | efetch -format xml > %s" % (i2q, wd+'data/'+tag+'_Efetch_Results_%s.xml' % (i)) for i2q,i in zip(chunks,range(len(chunks)))]
		xmldict=dict({})
		for i in range(len(chunks)):
			os.system(queries[i])
			with open(wd+'data/'+tag+'_Efetch_Results_%s.xml' % (i)) as xml:
				try:
					xmldict[i]=xmltodict.parse(xml.read())
				except Exception as e:
					print(e)
					print("there is a problem with one NCBI's XML file - discard batch %s" %(i))
		return xmldict		
		# [my_list[i * n:(i + 1) * n] for i in range((len(my_list) + n - 1) // n )]
		#hits=os.popen("cat test_out.txt | esearch -query '%s' | efetch -format uid" % (key)).read().split('\n')
	#return xmltodict.parse(hits[0])
#
xmlres=fetchlist(idlist,499)
###meshlist=[dpath.get(x,"**/MeshHeadingList/*") for x in xmlpmids]
xmlpmids=[dpath.search(xml,'**/PubmedArticle')['PubmedArticleSet']['PubmedArticle'] for xml in xmlres.values()] # (list of) list of dicts 
###
def parseml(xmlpmids):
	pmids=[dpath.get(x,"MedlineCitation/*PMID/*#text") for x in xmlpmids]
	pmids=[int(p) for p in pmids]
	### extract MESH
	meshlist=list(itertools.chain.from_iterable([dpath.values(x,"**/MeshHeadingList/*") for x in xmlpmids]))
	def meshdq(ref):
		descs=[dpath.values(m,"DescriptorName/**/#text") for m in ref]
		if len(descs)>0:
			quals=[dpath.values(m,"QualifierName/**/#text") for m in ref]
			# some descriptors have multiple qualifiers #
			quals=[q if len(q)>0 else ['None'] for q in quals]
			dqs=[["{}/{}".format(a_,b_) for a_, b_ in zip(d,q)] for d,q in zip(descs,quals)]
			dqs=list(set(itertools.chain.from_iterable(dqs)))
			dqs=[re.sub('/None','',s) for s in dqs]
			return dqs
		else: # no MeSH
			return []
	meshdict={p:meshdq(m) for p,m in zip(pmids,meshlist)}
	#	 
	### extract abstracts
	abslist=[dpath.values(x,"**/AbstractText") for x in xmlpmids]
	titlist=[dpath.values(x,"**/ArticleTitle") for x in xmlpmids] # [] if no title 
	#
	abslist=[a[0] if type(a)==list and len(a)>0 else a for a in abslist]
	#abslist=[list(set(itertools.chain.from_iterable(a))) for a in abslist]
	# "AbstractText" for those with no chapters, "#text" for the others	
	#abslist=[dpath.values(x,"**/AbstractText")+dpath.values(x,"**/#text") for x in abslist] # "AbstractText" for those with no chapters, "#text" for the others	
	def abschap(abst):
		# no abstract
		if type(abst)==list and len(abst)==0:
			return ''
		# AbstractText 
		elif type(abst)==str:
			return abst
		# list of ordered dicts #
		elif type(abst)==list:
			return ' '.join([dpath.get(od,'**/#text') for od in abst])
	#
	absform=[abschap(abst) for abst in abslist]
	# append title #
	absform=[' '.join(tit+[abst]) for tit,abst in zip(titlist,absform)]
	#absdict={p:abschap(abst) for p,abst in zip(pmids,abslist)}
	absdict={p:abst for p,abst in zip(pmids,absform)}
	return meshdict,absdict
##
if type(xmlpmids[0])==list:
	res_dicts=[parseml(xml) for xml in xmlpmids]
else:
	res_dicts=[parseml(xml) for xml in [xmlpmids]]
##
kw_dict=dict(ChainMap(*[t[0] for t in res_dicts]))
abs_dict=dict(ChainMap(*[t[1] for t in res_dicts]))
#######

abk=list(abs_dict.keys())

for k in abk:
	if abs_dict[k]=='':
		del abs_dict[k]


#
kwk=list(kw_dict.keys())

for k in kwk:
	if kw_dict[k]==[]:
		del kw_dict[k]		 


#
### STATISTICAL PRECAUTION: SET MUST BE IDENTICAL BETWEEN GENES AND MESH 
### I DON'T THINK THIS IS STILL NEEDED WITH THE CURRENT STATISTICS ###
# if "etype" not in globals():
# 	common=set(abs_dict.keys() & kw_dict.keys())# & set(tit_dict.keys()))
# elif not (etype.lower().startswith('g') or etype.lower().startswith('g')):
# 	common=set(abs_dict.keys() & kw_dict.keys())# & set(tit_dict.keys()))	
# else:
# 	common=set(abs_dict.keys()).union(set(kw_dict.keys()))
common=set(abs_dict.keys()).union(set(kw_dict.keys()))
# for d in [abs_dict,kw_dict,ful_dict,tit_dict]:
# 	dk=list(d.keys())
# 	for key in dk:
# 		if key not in common:
# 			del d[key]
# 		elif d[key]=='': # measure for abstracts
# 			del d[key]
# # recompute common
print("Effective set of PMIDs has a size of",len(common))
# get kws list
#abs_dict={k:v[:-1].split('!') for k,v in abs_dict.items()}
#kw_dict={k:v[:-1].split('!') for k,v in kw_dict.items()}
### REMOVE MESH TERMS THAT ARE PARENTS OF ANOTHER TERM IN THE SAME PAPER ### 
mesh_dict=decompress_pickle(sd+'input/meshD_names_to_ids.pbz2')
ful_dict=abs_dict
###
print("remove most-general, parent MeSHs for each PMID...")
def unparent(v):
	# remove MeSHs which aren't in MeSH dictionary - would have been removed anyway #
	v=[k for k in v if k.split('/')[0] in mesh_dict.keys()]
	# get a list of refence MeSH ids for the paper's keywords
	vsplit=[mesh.split('/')[0] for mesh in v]	
	vref=list(itertools.chain.from_iterable([mesh_dict[k] for k in vsplit]))
	v=[k for k in v if not any((itertools.chain.from_iterable([[ref.startswith(idk) for idk in mesh_dict[k.split('/')[0]]] for ref in set(vref)-set(mesh_dict[k.split('/')[0]])])))]
	return v
	
kw_dict={k:unparent(v) for k,v in kw_dict.items()} 
###
kws=list(set(itertools.chain.from_iterable(list(kw_dict.values()))))
print("number of keywords: %d" % len(kws))
#
print("create subsets of meshs")
#
keys_pub={k:[] for k in kws}
#
def shuffledict(old_it,new_dict=keys_pub):
	k,v=old_it
	for kw in v:
		new_dict[kw]+=[k]
	return k
#
pmids=set(map(shuffledict,list(kw_dict.items())))
### POSSIBILITY TO EXCLUDE ABBREVIATIONS ###
### load schartz-hearst algo ###
exec(open(sd+"code/schwartz_hearst.py").read())
#
ful_dict=dict({k:v for k,v in ful_dict.items() if k in common})
#abs_dict=dict({k:v for k,v in abs_dict.items() if k in pmids})
# cure abstracts
greek_alphabet = {
	b'\xce\x91':	 'Alpha',
	b'\xce\x92':	 'Beta',
	b'\xce\x93':	 'Gamma',
	b'\xce\x94':	 'Delta',
	b'\xce\x95':	 'Epsilon',
	b'\xce\x96':	 'Zeta',
	b'\xce\x97':	 'Eta',
	b'\xce\x98':	 'Theta',
	b'\xce\x99':	 'Iota',
	b'\xce\x9a':	 'Kappa',
	b'\xce\x9b':	 'Lamda',
	b'\xce\x9c':	 'Mu',
	b'\xce\x9d':	 'Nu',
	b'\xce\x9e':	 'Xi',
	b'\xce\x9f':	 'Omicron',
	b'\xce\xa0':	 'Pi',
	b'\xce\xa1':	 'Rho',
	b'\xce\xa3':	 'Sigma',
	b'\xce\xa4':	 'Tau',
	b'\xce\xa5':	 'Upsilon',
	b'\xce\xa6':	 'Phi',
	b'\xce\xa7':	 'Chi',
	b'\xce\xa8':	 'Psi',
	b'\xce\xa9':	 'Omega',
	b'\xce\xb1':	 'alpha',
	b'\xce\xb2':	 'beta',
	b'\xce\xb3':	 'gamma',
	b'\xce\xb4':	 'delta',
	b'\xce\xb5':	 'epsilon',
	b'\xce\xb6':	 'zeta',
	b'\xce\xb7':	 'eta',
	b'\xce\xb8':	 'theta',
	b'\xce\xb9':	 'iota',
	b'\xce\xba':	 'kappa',
	b'\xce\xbb':	 'lamda',
	b'\xce\xbc':	 'mu',
	b'\xce\xbd':	 'nu',
	b'\xce\xbe':	 'xi',
	b'\xce\xbf':	 'omicron',
	b'\xcf\x80':	 'pi',
	b'\xcf\x81':	 'rho',
	b'\xcf\x83':	 'sigma',
	b'\xcf\x84':	 'tau',
	b'\xcf\x85':	 'upsilon',
	b'\xcf\x86':	 'phi',
	b'\xcf\x87':	 'chi',
	b'\xcf\x88':	 'psi',
	b'\xcf\x89':	 'omega',}
#
def greek_sub_a(m):
	if m.group() is not None:
		try: 
			return greek_alphabet[m.group().encode('utf-8')]
		except:
			return ''
#
def greek_sub_f(m):
	if m.group() is not None:
		try: 
			return greek_alphabet[m.group().encode('utf-8')]
		except:
			return ' '
#
def abb_pairs(s):
	random.seed(2202)
	d=extract_abbreviation_definition_pairs(doc_text=s,first_definition=True)
	d={str(k):str(v) for k,v in d.items()}
	# Deal with abbreviations of plurals etc...
	pattern = re.compile("^[a-z][A-Z]+$|^[A-Z]+[a-z]$")
	items=list(d.items())
	for k,v in items:
		if pattern.match(k):
			#print("found plural/unlemmatized abbreviation:",k)
			sing=re.sub("^[a-z]|[a-z]$","",k)
			#print("lemmatization:",sing)
			d[sing]=v
	return d #{str(k):str(v) for k,v in d.items()}
#
def celliNER(items,model):
	random.seed(2202)
	nlp = globals()[model]
	items = [(x[1],x[0]) for x in items]
	spacy.util.fix_random_seed(2202)
	docs = nlp.pipe(items, as_tuples=True)
	d={int(tpl[1]):list(set([str(ent) for ent in tpl[0].ents if ent.label_ in ["CELL_TYPE", "CELL_LINE"]])) for tpl in docs}
	return d
#
def distro(data,ncores,d1):
	chunks=[x.tolist() for x in np.array_split(data,ncores)]
	d1s=[d1+str(i) for i in range(0,ncores)]
	return iter([(chunks[i],d1s[i]) for i in range(0,ncores)])
#
# find cell line names 
print("loading scispacy model")
nlp = spacy.load("en_ner_jnlpba_md")

if len(list(ful_dict.items())) < ncores:
	ncores = len(list(ful_dict.items()))

var_names_nlp=["nlp_"+str(i) for i in range(ncores)]
for name in var_names_nlp:
	globals()[name] = nlp

with Pool(ncores) as pool:
	print("scispacy starmap ...")
	res=pool.starmap(celliNER,distro(list(ful_dict.items()),ncores,'nlp_'))

cell_dict={}
for d in res:
	cell_dict.update(d)


# cell_dict = {key:nlp(val).ents for key,val in ful_dict.items()}
# cell_dict = {key:list(set([str(ent) for ent in cell_dict[key] if ent.label_ in ["CELL_TYPE", "CELL_LINE"]])) for key in list(cell_dict.keys())}


#abs_dict=dict({k:[re.sub('[\u0080-\uFFFF]',greek_sub_a,i) for i in v] for k,v in abs_dict.items()})
print("finding abbreviations")
ful_dict=dict({k:re.sub('[\u0080-\uFFFF]',greek_sub_f,v) for k,v in ful_dict.items()})
abb_dict=dict({key:abb_pairs(val) for key,val in ful_dict.items()})
abb_dict={pmid:{ab:df for ab,df in abb_dict[pmid].items() if df not in cell_dict[pmid]} for pmid in list(abb_dict.keys())} 
#abb_dict=pd.Series(abb_dict)
#print(ful_dict)
#print(abb_dict)
keys_pub={k:list(set(v)) for k,v in keys_pub.items()}
#
# print("Construct input file for D3NER")
# #
# def d3ner(i,tit_dict=tit_dict,ful_dict=ful_dict):
# 	if tit_dict[i]!='' and ful_dict[i]!='':
# 		S=str(i)+'|t|'+tit_dict[i]+'\n'+str(i)+'|a|'+ful_dict[i]+'\n\n'
# 		return S
# #
# d3=[d3ner(i) for i in ful_dict.keys()]
# d3=[i for i in d3 if i]
# #
# S=''.join(d3)
#
#print("Merge titles and abstracts")
#
#ful_dict={i:tit_dict[i]+ful_dict[i] for i in ful_dict.keys()}
#
print("save dictionaries")
#
#compress_pickle(wd+'data/'+tag+'_ids_to_abs', abs_dict)
### take into account the entity type selection ### 
if 'etype' in globals():
	if etype.lower().startswith('m'):
		### then no gene extraction ###
		compress_pickle(wd+'data/'+tag+'_ids_to_kw', keys_pub)
		compress_pickle(wd+'data/'+tag+'_ids_to_abbr', dict({}))
		compress_pickle(wd+'data/'+tag+'_ids_to_cells', dict({}))
		compress_pickle(wd+'data/'+tag+'_ids_to_full', dict({}))
	elif etype.lower().startswith('g'):		
		### then no MeSH extraction ###
		compress_pickle(wd+'data/'+tag+'_ids_to_kw', dict({}))
		compress_pickle(wd+'data/'+tag+'_ids_to_abbr', abb_dict)
		compress_pickle(wd+'data/'+tag+'_ids_to_cells', cell_dict)
		compress_pickle(wd+'data/'+tag+'_ids_to_full', ful_dict)
	else:
		print("### WARNING: unrecognised etype: %s ###\n### did You specify your choice (gene/MeSH) correctly? ###\n### proceed to ignore etype ###" % (etype))		
		compress_pickle(wd+'data/'+tag+'_ids_to_kw', keys_pub)
		compress_pickle(wd+'data/'+tag+'_ids_to_abbr', abb_dict)
		compress_pickle(wd+'data/'+tag+'_ids_to_cells', cell_dict)
		compress_pickle(wd+'data/'+tag+'_ids_to_full', ful_dict)
else:
	compress_pickle(wd+'data/'+tag+'_ids_to_kw', keys_pub)
	compress_pickle(wd+'data/'+tag+'_ids_to_abbr', abb_dict)
	compress_pickle(wd+'data/'+tag+'_ids_to_cells', cell_dict)
	compress_pickle(wd+'data/'+tag+'_ids_to_full', ful_dict)
###
print("done")