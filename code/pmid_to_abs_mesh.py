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
from multiprocess import Process, Manager
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
#
english_words_set= english_words_set.union(set(words.words()))
#
Entrez.email='Luca.musella@uk-erlangen.de'
Entrez.api_key='b89451c47164b8705c76108eec6386658607'
#
wd= sys.argv[1] #'/home/musellla/tam_textmining/'
tag= sys.argv[2]
idlist= sys.argv[3]
#
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
	b'\xcf\x89':	 'omega',
}

# special handling of greek letters: remove them from english dictionary

english_words_set=english_words_set - set(greek_alphabet.values())	

#
#load batch of ids (a .txt file)
# decompress_pickle(wd+'data/pmid_queries_to_efetch.pbz2')
#idlist=idlist[:20]
with open(idlist,"r") as file:
	idlist=file.read().splitlines() #file.readlines()
#
idlist=[int(i) for i in idlist]
#
print("list of PMIDs loaded")

id_abs={}


def chiave(pmid,lab,d):
	#print(pmid,', ',lab)
	x=read(efetch(db='pubmed', id=pmid, retmode='xml', apikey='b89451c47164b8705c76108eec6386658607'))
	x=x['PubmedArticle'][0]
	xdf=pd.DataFrame(x)
	article=xdf.loc[['Article'],['MedlineCitation']].values[0]
	article=article.tolist()[0]
	try:
		abs_text=article['Abstract']['AbstractText']
		if len(abs_text)>1:
			abs_text=[str(s) for s in abs_text]
			abs_text=' '.join(abs_text)
		else:
			abs_text=abs_text[0]
	except:
		print("No abs?",pmid)
		abs_text=' '		
	#print(pmid,', ',lab)	
	#creating dict with pmid:abs
	d[pmid]=abs_text
	#print(pmid)#
	#print(abs_text)
	# you check the presence of keywords by
	# u'KeywordList' in x['MedlineCitation'] OR
	# u'MeshHeadingList' in x['MedlineCitation']
	if u'KeywordList' in x['MedlineCitation'] or u'MeshHeadingList' in x['MedlineCitation']:
		try:
			kws=list(xdf.loc['KeywordList']['MedlineCitation']) # it's a "list element"
		except:
			kws=[]	
		try: # this is the optimal case
			mesh=xdf.loc['MeshHeadingList']['MedlineCitation'] # it's a list #
			mesh_q=[list(q['QualifierName']) for q in mesh]
			mesh_d=[str(q['DescriptorName']) for q in mesh]
			mesh_q=[[q.title() for q in Q] for Q in mesh_q] 
			[l.append(d) for l,d in zip(mesh_q, mesh_d)]
			mesh_q=[[q.lower() for q in Q] for Q in mesh_q] 
			mesh_qt=[list(itertools.combinations(q,min(len(q),2))) for q in mesh_q]
			mesh_qt=list(itertools.chain.from_iterable(mesh_qt))
			mesh_qt=list(set(mesh_qt))
			#print(mesh_qt)
			mesh_qt=[(pmid,lab,'MESH')+q for q in mesh_qt]
			if len(kws) > 0:
				kws=[str(i) for i in list(kws[0])] 
				meshes=list(itertools.chain.from_iterable(mesh_q))
				kws=[i.lower() for i in kws if i not in meshes]
				kws=list(itertools.combinations(kws,2))
				kws=[(pmid,lab,'KW')+q for q in kws]
				mesh_qt=mesh_qt+kws
			return mesh_qt
#merge keywords and mesh words
		except: # this runs if the article has no mesh terms
			if len(kws) == 0:
				return [(pmid,lab,'MESH/KW','NA','NA')]
			else: # this runs if there is no mesh but some general keywords
				mesh_qt=[(pmid,lab,'MESH','NA','NA')]
				kws=[str(i).lower() for i in list(kws[0])]
				kws=list(itertools.combinations(kws,2))
				kws=[(pmid,lab,'KW')+q for q in kws]
				mesh_qt=mesh_qt+kws
				return mesh_qt
	else: # nothing to retrieve anyway
		return [(pmid,lab,'MESH/KW','NA','NA')]

#ip(tam_id,tam_tx)}
print("fetch PMIDs")
i=0
pmid_to_kw=[]
while i<len(idlist):
	#print("counter: ",i)
	try:
		res=chiave(idlist[i],"tam[mesh]",id_abs)
		pmid_to_kw.append(res)
		i+=1
		print("efetch: %.2f percent done" % (100*(i/len(idlist))))
		#time.sleep(3)
	except Exception as e:
		print("Error! Bad internet connection?")
		print(e)
		time.sleep(5)
#

# compress_pickle(wd+'data/iter1_ids_to_kw', pmid_to_kw)
# compress_pickle(wd+'data/iter1_ids_to_abs', id_abs)

#################
################# code to use for converting pmid_to_kw format

# # load efectch results 
# #
# ### tam[mesh] papers only
# id2kw=decompress_pickle(wd+"data/"+tag+"_ids_to_kws.pbz2")
####
# process and save IDs to MESH dictionaries

pmid_to_kw=list(itertools.chain.from_iterable(pmid_to_kw))
df_TP=pd.DataFrame(pmid_to_kw,columns=['PMID','LABEL','SOURCE','KW_1','KW_2'])
#
df=df_TP.drop_duplicates(ignore_index=True)
# RESTRICT TO MESH TERMS ONLY
df=df[df['SOURCE']=='MESH']
#

#
del df_TP
#
pmids=df['PMID']
pmids=list(set(pmids))
# moreover, reduce to mesh terms
test=df[df['PMID'].isin(pmids)]
print("data loaded and preprocessed")
# WE NEED:
# 1) PMID-TO-KWS MAP
# 2) LIST OF ALL KWS IN EVERY PAPER
# I will obtain a matrix o 0s and 1s, PMIDxKW
# Compute distance matrix 
# Hierarchical clustering
# "Cut" the tree at a given density/leve;
#
def pmid_to_kws(pmid,df):
	df=df[df['PMID']==pmid]
	lab=df[['LABEL']].iloc[0].values[0]
	kw1=list(df['KW_1'])
	kw2=list(df['KW_2'])
	kws=list(set(kw1+kw2))
	kws=[i for i in kws if i]
	return (pmid,lab,kws)
kws_in_pmid=[pmid_to_kws(pm,test) for pm in pmids]
kws=[t[2] for t in kws_in_pmid]
kws=list(set(list(itertools.chain.from_iterable(kws))))
kws.sort()
#
print("number of keywords: %d" % len(kws))
#
print("create subsets of meshs")
# dictionarize kws_to_pubmed
keys_pub=dict({key:[] for key in kws})
for key in keys_pub:
	keys_pub[key].append([t[0] for t in kws_in_pmid if key in t[2]])
	keys_pub[key]=list(set(itertools.chain.from_iterable(keys_pub[key])))
#
print("save MeSH:PMIDs dictionary")
#
### STATISTICAL PRECAUTION: SET MUST BE IDENTICAL BETWEEN GENES AND MESH 
common=list(keys_pub.keys() & id_abs.keys())
for d in [keys_pub,id_abs]:
	intd=dict({k:d[k] for k in common})
	d=intd
###
#compress_pickle(wd+'data/'+tag+'_ids_to_kw_unprocessed', pmid_to_kw)
compress_pickle(wd+'data/'+tag+'_ids_to_kw', keys_pub)
# NOT WORKING compress_pickle(wd+'data/'+tag+'_ids_to_abs_raw', id_abs)
# ID_TO_ABS_RAW AS .JSON.gz
# remove greek letters
def greek_sub(m):
	if m.group() is not None:
		try: 
			return greek_alphabet[m.group().encode('utf-8')]
		except:
			return None

id_abs={k:re.sub('[\u0080-\uFFFF]',greek_sub,v) for k,v in id_abs.items()}
#
with gzip.open(wd+'ids_to_abs_raw.json.gz', 'at', encoding='UTF-8') as zipfile:
    json.dump(id_abs, zipfile)
#
#############
#############
# process and save IDs to Abs dictionaries
def cleansing_1x(x,english_words_set=english_words_set):
	# Works on one abs at a time
	x=x.replace('\n',' ')
	x=re.sub('\/',' ',x)
	x=re.sub(',',' ',x)
	x=x.replace('.','')
	x=re.sub('[^\\w-]+',' ',x)
	# def re_dash(m):
	# 	# there is at most one dash per gene (?)
	# 	def window(seq, n=2):
	# 	#     "Returns a sliding window (of width n) over data from the iterable"
	# 	# 		"   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
	# 		it = iter(seq)
	# 		result = tuple(islice(it, n))
	# 		if len(result) == n:
	# 			yield result    
	# 		for elem in it:
	# 			result = result[1:] + (elem,)
	# 			yield result
	# #
	# 	def allsplits(strng, sep):
	# 		strng = strng.split(sep)
	# 		return [sep.join(strng[:i]) for i in range(1,len(strng)+1)]
	# 	if m.group() is not None:
	# 		s1= window(re.split('\-',m.group()),n=1) #allsplits(m.group(),'-') #re.split('\-',m.group())#re.split('\-',m.group()) #
	# 		s1=["".join(t) for t in s1]			
	# 		s1=' '.join(s1)
	# 		s2=window(re.split('\-',m.group()))
	# 		s2=["-".join(t) for t in s2]			
	# 		s2=' '.join(s2)
	# 		return s1+s2 # m.group()+s
	# x=re.sub(" [A-Za-z\d\u0080-\uFFFF]*\-[A-Za-z\d\u0080-\uFFFF-]* ", re_dash,x)
	#
	# find all non-ascii, manipulate in separate list and concatenate
	def find_greek(ab,greek_alphabet=greek_alphabet):
		def greek_sub(m):
			if m.group() is not None:
				try: 
					return greek_alphabet[m.group().encode('utf-8')]
				except:
					return None
		def undash(w):
			notdash=''.join(w.split('-'))
			return w+' '+notdash
		mat=re.findall('[a-zA-Z\d\-]+[\u0080-\uFFFF][a-zA-Z\d]*|[\u0080-\uFFFF]+\-[a-zA-Z\d]+',ab)
		if len(mat)==0:
			return ab
		else:
			# version without any greek letter -> WRONG, produces artifacts! (IL-1beta->IL-1)
			
			#matout=[re.sub('[\u0080-\uFFFF]','',i) for i in mat]
			#matout=' '.join(matout)
			#
			mat=list(set(mat))
			# translation of greek characters to strings
			mat=[re.sub('[\u0080-\uFFFF]',greek_sub,i) for i in mat]
			#
			mat=[undash(w) for w in mat]
			mat=' '.join(mat)
			return ab+' '+mat+' '#+matout
	#
	x=find_greek(x)
	# refinements
	x=re.sub(' [\u0080-\uFFFF-\d]* ', ' ', x)
	x=re.sub(r' \d+ ','', x)
	# multiple spaces become single spaces
	x=re.sub(r'\s+', ' ', x)
	# LEMMATIZATION
	stemmer = WordNetLemmatizer()
	def lemma(document,mod):
		if mod=='a':
			document=document.split()
			document=[stemmer.lemmatize(w, pos='a') for w in document]
			document=' '.join(document)
		if mod=='v':
			document=document.split()
			document=[stemmer.lemmatize(w, pos='v') for w in document]
			document=' '.join(document)
		else:
			document=document.split()
			document=[stemmer.lemmatize(w) for w in document]
			document=' '.join(document)
		return document	
	#
	x=lemma(x,'n')
	x=lemma(x,'v')
	x=lemma(x,'a')
	#	
	# REMOVE ALL ENGLISH WORDS 
	x=x.split()
	x=[w for w in x if w.lower() not in english_words_set]
	#
	def clean_dash(w,P,english_words_set=english_words_set):
		def window(seq, n=2):
		# "Returns a sliding window (of width n) over data from the iterable"
		# s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
			it = iter(seq)
			result = tuple(islice(it, n))
			if len(result) == n:
				yield result    
			for elem in it:
				result = result[1:] + (elem,)
				yield result
		# pattern "-word" or "word-"
		ds=list(window(w.split('-')))
		if len(ds)==0:
			return w
		else:
			# no "word-word-word" is possible
			for pair in ds:
				pair=[stemmer.lemmatize(p, pos='a') for p in pair]
				pair=[stemmer.lemmatize(p, pos='v') for p in pair]
				pair=[stemmer.lemmatize(p, pos='n') for p in pair]
				if not all(p.lower() in english_words_set for p in pair):
					# then remove the english word
					g='-'.join(pair)
					P.append(g)	
	P=[]
	x=[clean_dash(w,P) for w in x]	
	x+=P
	x=[w for w in x if w]
	x=[w for w in x if w != '']		
	# prepare stopwords
	all_stopwords = stopwords.words('english')
	all_stopwords+=[i.upper() for i in all_stopwords]
	all_stopwords+=[i.capitalize() for i in all_stopwords]
	all_stopwords+=list(string.ascii_lowercase)
	all_stopwords+=list(string.ascii_uppercase)
	# stop "*co*-occurrence"
	all_stopwords+=['co','Co','CO','M1','M2']
	#
	x=list(set(x))
	# more careful evaluation with lowercase (not necessary if you use reference aliases)
	# x = [i.lower() for i in x]
	return '!'.join(x)

print("Cleaning abstracts...")
test={}
for key in list(id_abs.keys()):
	test[key]=cleansing_1x(id_abs[key])	


# set and sort, stopwords
 # macrophage states are buzzwords
#
def set_and_sort(txt):
	return ' '.join(set(txt.split()))

#
for key in test:
	test[key]=set_and_sort(test[key])		

#
def stop(doc,all_stopwords):
	doc=doc.split()
	return [i for i in doc if not i in all_stopwords]

#
for key in test:
	test[key]=stop(test[key],all_stopwords)		

#
pmid_to_abs=test

compress_pickle(wd+'data/'+tag+'_ids_to_abs', pmid_to_abs)
#

#test={key:cleansing(val) for key,val in list(id_abs.items())[:20]}
