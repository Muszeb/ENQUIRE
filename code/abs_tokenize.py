import os
import sys
import numpy as np 
import pandas as pd
import itertools
from itertools import islice
import multiprocessing as mp
from multiprocessing import Pool
from collections import OrderedDict
from multiprocess import Process, Manager
import re
import random
import time
from datetime import datetime
import bz2
import json
import threading
import string
from english_words import english_words_set
from Bio import pairwise2
from itertools import starmap
import csv
import subprocess
import nltk
nltk.download('stopwords')
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
nltk.download('words')
from nltk.corpus import words
from nltk.stem import WordNetLemmatizer
import spacy
from scispacy.abbreviation import AbbreviationDetector
try:
    import _pickle as pickle
except ModuleNotFoundError:
    	import pickle

# Pickle a file and then compress it into a file with extension 
def compress_pickle(title, data):
	with bz2.BZ2File(title + ".pbz2", "w") as f: 
		pickle.dump(data, f)


# Pickle a file and then compress it into a file with extension 
def decompress_pickle(file):
	data = bz2.BZ2File(file, "rb") 
	data = pickle.load(data)
	return data

#
wd= sys.argv[1] # tmp
sd= re.split("tmp",wd)[0]
# wd='/home/musellla/tam_textmining/'
tag= sys.argv[2] #iter1
#
abb_dict=decompress_pickle(wd+'/data/'+tag+'_ids_to_abbr.pbz2')
ful_dict=decompress_pickle(wd+'data/'+tag+'_ids_to_full.pbz2')
cell_dict=decompress_pickle(wd+'data/'+tag+'_ids_to_cells.pbz2')
#
print("modules and data loaded")
#
if all([len(d)==0 for d in [abb_dict,ful_dict,cell_dict]]):
	print("### WARNING: all necessary dictionary have length 0 - did you set etype=Mesh? ###")
	compress_pickle(wd+'data/'+tag+'_ids_to_abs', dict({}))
else:
	print("remove abbreviation definitions from text")
	#
	def absurge(i,ful=ful_dict,abb=abb_dict,cell=cell_dict):
		#0.
		if abb[i] or cell[i]:
			# LIST ORDER MAY DIFFER BETWEEN RUNS! DIFFERENT RESULTS COULD BE OBTAINED WHEN REMOVING STRINGS #
			# use alphabetica, THEN decreasing length (reverse) as a proxy for more specific -> less specific order
			for v in sorted(sorted(list(abb[i].values()) + cell[i]),key=len,reverse=True):
				ful[i]=re.sub(re.escape(v), "", ful[i])
	#
	for i in ful_dict.keys():
		absurge(i)
	#
	### TOKENIZATION ### 
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
	english_words_set= english_words_set.union(set(words.words()))
	english_words_set=english_words_set - set(greek_alphabet.values())	
	#
	#############
	print("tokenize...")
	# process and save IDs to Abs dictionaries
	def cleansing_1x(i,ful,english_words_set,greek_alphabet):
		#
		ful=globals()[ful]
		english_words_set=globals()[english_words_set]
		greek_alphabet=globals()[greek_alphabet]
		# Works on one abs at a time
		x=ful[i]
		#
		x=x.replace('\n',' ')
		x=re.sub('\/',' ',x)
		x=re.sub(',',' ',x)
		x=x.replace('.',' ')
		x=re.sub('[^\\w-]+',' ',x)
		x=re.sub(r' \d+ ',' ', x)
		# multiple spaces become single spaces
		x=re.sub(r'\s+', ' ', x)
		# find all non-ascii, manipulate in separate list and concatenate
		def find_greek(ab,greek_alphabet=greek_alphabet):
			def greek_sub(m):
				if m.group() is not None:
					try: 
						return greek_alphabet[m.group().encode('utf-8')]
					except:
						return ' '
			def undash(w):
				notdash=''.join(w.split('-'))
				return w+' '+notdash
			mat=re.findall('[a-zA-Z\d\-]*[\u0080-\uFFFF][a-zA-Z\d]*|[\u0080-\uFFFF]+\-[a-zA-Z\d]+',ab)
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
		x=re.sub('\*', ' ', x)
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
		#
		def clean_dash(w,P,english_words_set=english_words_set):
			
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
		#
		def dedash(w,P,english_words_set=english_words_set):
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
			def lemma(pair):
				pair=[stemmer.lemmatize(p, pos='a') for p in pair]
				pair=[stemmer.lemmatize(p, pos='v') for p in pair]
				pair=[stemmer.lemmatize(p, pos='n') for p in pair]
				return(pair)
			#
			ws=w.split('-')
			if len(ws)==1:
				return w
			else:
				ws=lemma(ws)
				if ws[-1] in english_words_set: # "gene-related", "gene-dependent". no .lower(), because always comes after at least 1 word!
						del ws[-1]
				if ws[0].lower() in english_words_set and len(ws[0])>2: # "pre-caspase", "anti-PDCD1"
						del ws[0]
				if len(w.split('-'))>2:
					ds=list(window(w.split('-')))
					ds=[t for t in ds if not all([j in english_words_set for j in t])]
					ds=['-'.join(list(t)) for t in ds]
					P+=ds
				return '-'.join(ws)
		P=[]
		x=[dedash(w,P) for w in x]	
		x+=P
		x=[w for w in x if w.lower() not in english_words_set]
		x=[w for w in x if w]
		x=[w for w in x if w != '']		
		x=[re.sub("^-|-$",'',w) for w in x]
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
		return i,x #'!'.join(x)
	#
	def greek_sub_a(m):
		if m.group() is not None:
			try: 
				return greek_alphabet[m.group().encode('utf-8')]
			except:
				return ''
	#
	print("duplicate dictionaries for parallelization...")

	ncores= int(sys.argv[3]) #32

	var_names_abs=["ful_dict_"+str(i) for i in range(ncores)]
	var_names_syn=["english_words_set_"+str(i) for i in range(ncores)]
	var_names_abb=["greek_alphabet_"+str(i) for i in range(ncores)]

	for name in var_names_abs:
		globals()[name] = ful_dict	


	for name in var_names_syn:
		globals()[name] = english_words_set #aliasd=pd.Series(aliasdf['aliases'].tolist(),index=aliasdf.symbol).to_dict()


	for name in var_names_abb:
		globals()[name] = greek_alphabet
	#
	def distro(datal,ncores,d1,d2,d3):
		# d1 and d2 are strings
		# set repetitions
		rep=int(len(datal)/ncores)+2
		d1s=[d1+str(i) for i in range(0,ncores)]*rep
		d2s=[d2+str(i) for i in range(0,ncores)]*rep
		d3s=[d3+str(i) for i in range(0,ncores)]*rep
		#print(d1s)
		#d2s=[d2+str(i) for i range(0,ncores)]*rep
		return iter((datal[i],d1s[i],d2s[i],d3s[i]) for i in range(0,len(datal)))
	#
	with Pool(ncores) as pool:
		#res=[pool.apply(find_in_series,args=(idd,globals()["abs_dict_"+str(ids.index(idd))],globals()["syn_dict_"+str(ids.index(idd))])) for idd in ids]   
		print("starmap...")
		abs_dict=dict(pool.starmap(cleansing_1x,distro(list(ful_dict.keys()),ncores,'ful_dict_','english_words_set_','greek_alphabet_')))
		print("store in dictionary")
		#res=dict(res)
	#
	abs_dict=dict({k:[re.sub('[\u0080-\uFFFF]',greek_sub_a,i) for i in v] for k,v in abs_dict.items()})
	#
	compress_pickle(wd+'data/'+tag+'_ids_to_abs', abs_dict)