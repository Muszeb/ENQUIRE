# gene extraction by writing to file (maybe achieves perfect parallelization) 
# exec(open("/home/musellla/tam_textmining/code/gene_extraction.py").read())
import os
import sys
from nltk.corpus.reader.plaintext import PlaintextCorpusReader
import numpy as np 
import pandas as pd
import itertools
from itertools import islice
from collections import OrderedDict
import re
import random
import time
from datetime import datetime
import nltk
nltk.download('stopwords')
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
import networkx as nx
import bz2
import json
import threading
from multiprocess import Process, Manager
import multiprocessing as mp
from multiprocessing import Pool
import string
from scipy.spatial import distance_matrix
from Bio import pairwise2
nltk.download('words')
from nltk.corpus import words
from nltk.stem import WordNetLemmatizer
from english_words import english_words_set
from Bio import pairwise2
from itertools import starmap
import csv
import warnings
try:
    import _pickle as pickle
except ModuleNotFoundError:
    	import pickle

#
warnings.simplefilter(action='ignore', category=FutureWarning)
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
def opensplitFile(path):
	#param path: path/to/file.ext (str)
	#Returns contents of file (str)
	with open(path) as file:
		data = file.read().split('\n\n\n')
	return data
#
#sys.argv=['', '/home/musellla/txtmining/tmp-MEL_D3NER/MEL_D3NER_subgraphs_expansion1/', 'MEL_D3NER_subgraphs_expansion1', '32']
#
wd= sys.argv[1] # tmp
sd= re.split("tmp",wd)[0]
# wd='/home/musellla/tam_textmining/'
tag= sys.argv[2] #iter1
#
### load data ###
pmid_to_abs = decompress_pickle(wd+'/data/'+tag+'_ids_to_abs.pbz2')
abb_dict = decompress_pickle(wd+'/data/'+tag+'_ids_to_abbr.pbz2')
###
if len(pmid_to_abs)==0:
	print("PMID to abstract dictionary is empty - skip NER")
	compress_pickle(wd+"data/"+tag+"_ids_to_genes",dict({}))
###
else:
	english_words_set= english_words_set.union(set(words.words()))
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
	# SETUP GENE DICTIONARIERS FOR H.s and M.m
	# human and mouse genes
	#aliasdf = pd.read_csv(sd+'input/gene_aliases_df.tsv',sep='\t')
	aliasdf = pd.read_csv(sd+'input/aliases_uniprotensembl_as_df_Paralogues.tsv',sep='\t')
	aliasdf=aliasdf.dropna()
	#
	def alclean(al):
		def greek_sub_a(m):
			if m.group() is not None:
				try: 
					return greek_alphabet[m.group().encode('utf-8')]
				except:
					return ''
			#
		if al in english_words_set:
			return False
		elif bool(re.search("^[a-zA-Z]{1,2}$|^[\W\d]+$", al)):
			return False
		elif al in greek_alphabet.values():
			return False
		else:
			al=re.sub("[\(\)\[\]]",'',al)
			al=re.sub('[\u0080-\uFFFF]',greek_sub_a,al)
			return al
	#
	aliasdf['alias_symbol']=aliasdf['alias_symbol'].apply(alclean,1)
	# ``ser.astype(object).apply()``
	# remove 1-2 letters only aliases
	aliasdf=aliasdf[~(aliasdf.alias_symbol==False)]
	# aliases=aliasdf.aliases.tolist()
	# # unless they are reference genes => regardless being a reference gene
	# #refs=aliasdf.symbol.tolist()
	# al=al.split(';')
	# al=[a for a in al if a not in english_words_set]
	# al=[a for a in al if not bool(re.search("^[a-zA-Z]{1,2}$|^[\W\d]+$", a))]
	# al=[re.sub("[\(\)\[\]]",'',a) for a in al]
	# al=[a for a in al if a not in greek_alphabet.values()]	
	# return ';'.join(al)
	# # #
	# aliasdf['aliases']=aliasdf.aliases.apply(alclean,1)
	# #
	# # convert to an alias-symbol series
	# def serialal(aliasdf):
	# 	def reser(row):
	# 		als=row[1].split(';')
	# 		sym=row[0]
	# 		return list(itertools.product(als,[sym]))	
	# 	#
	# 	xl=aliasdf.apply(reser,1)
	# 	nl=list(itertools.chain(*xl.tolist()))
	# 	nl=pd.DataFrame(nl)
	# 	nl.index=nl[0]
	# 	AliasS=pd.Series(nl[1])
	# 	return AliasS
	# #
	#aliass=serialal(aliasdf)	
	aliass=pd.Series(aliasdf.symbol.tolist(),index=aliasdf.alias_symbol.tolist())
	#aliass=aliass[~aliass.index.duplicated(keep='first')]
	#
	#dictionarize
	#
	#bad_hits=list(set(filter(lambda x: re.match("^[a-zA-Z]{1,2}$",str(x)),aliases)))
	#aliasdf=aliasdf[~aliasdf['aliases'].isin(bad_hits)]
	#
	syms=aliasdf.symbol.tolist()
	###
	# with gzip.open(wd+'ids_to_abs_raw.json.gz', 'rt', encoding='UTF-8') as zipfile:
	#     pmids_to_raw_abs = json.load(zipfile)
	# #
	# # ACHTUNG! KEYS ARE NOW STRINGS
	# pmids_to_raw_abs={int(k):v for k,v in pmids_to_raw_abs.items()} 
	# # sub greek letters

	###
	print("PMIDs and alias DF loaded")
	#
	print("set multithreading")
	##
	print("duplicate PMIDs-abstracts and aliasases dictionary for parallelization...")

	ncores= int(sys.argv[3]) #32

	var_names_abs=["abs_dict_"+str(i) for i in range(ncores)]
	var_names_syn=["syn_dict_"+str(i) for i in range(ncores)]
	var_names_abb=["abb_dict_"+str(i) for i in range(ncores)]

	for name in var_names_abs:
		globals()[name] = pmid_to_abs	


	for name in var_names_syn:
		globals()[name] = aliass #aliasd=pd.Series(aliasdf['aliases'].tolist(),index=aliasdf.symbol).to_dict()


	for name in var_names_abb:
		globals()[name] = abb_dict
	#
	# to avoid overtaking the max number of threads, split the list (feed as first "climbME" argument)
	genes_pub=dict({sym:[] for sym in syms})
	#chunks = iter([list(pmid_to_abs.keys())[x:x+ncores] for x in range(0, len(list(pmid_to_abs.keys())), ncores)])
	#
	# bingo=['OAS1',
	# 	'Skull21',
	# 	'Skull14',
	# 	'Skull22',
	# 	'IGKV1-16',
	# 	'Sdtq10',
	# 	'Rpl18',
	# 	'Rpl14',
	# 	'Rpl13',
	# 	'CYRIB',
	# 	'Skull20',
	# 	'Skull16',
	# 	'Lgals1',
	# 	'Sdtq11',
	# 	'L1cam',
	# 	'Pdilt',
	# 	'PDCD1']
	#	


	# def find_in_tokens(idd,abs_dict,syn_dict):
	# 	#
	# 	tok=abs_dict[idd]
	# 	#
	# 	def mappal(alit):
	# 	#print(alit)
	# 		k,v = alit
	# 		#bool(re.search("[\W\s]+al+[\W\s]",txt))
	# 		#print(v)
	# 		if any(al in tok for al in v.split(';')):
	# 			#print("found",k)
	# 			return k
	# 	#
	# 	hit=list(map(mappal,syn_dict.items())) 
	# 	# df[df.columns.intersection(l)]
	# 	hit=[i for i in hit if i]
	# 	# hit contains gene entities. What about "GENE1-GENE2" entities (if any)?
	# 	# delete any word that matches an alias per se (if so, should be in hit already)
	# 	allal=';'.join(syn_dict.values()).split(';')
	# 	#
	# 	g1g2= list(set(tok) - set(allal))
	# 	g1g2=[w for w in g1g2 if '-' in w]		
	# 	#
	# 	if len(g1g2) > 0:
	# 		#	
	# 		#print("A gene-gene occurrence(?)",g1g2)
	# 		for pair in g1g2:
	# 			pair=pair.split('-')
	# 			# BOTH must be gene aliases ("GENE-dirt" occurrences will be discarded)
	# 			if all(p in allal for p in pair):
	# 				#
	# 				for p in pair:
	# 				#print("found", pair[0],'-',pair[1],'pair...')
	# 				# assuming an alias is unique to its reference gene:
	# 					hit+=[k for k in syn_dict if p in syn_dict[k].split(';')] 
	# 	hit=list(set(hit))
	# 	hit=[(h,idd) for h in hit]
	# 	return hit 

	# one chunk multithreaded at a time

	#
	def find_in_series(idd,abs_dict,syn_dict,abb_dict):
		#
		#print(idd,abs_dict,syn_dict,abb_dict)
		abs_dict=globals()[abs_dict]
		syn_dict=globals()[syn_dict]
		abb_dict=globals()[abb_dict]
		#
		tok=abs_dict[idd]
		abb=abb_dict[idd]
		#print(tok,abb)
		#
		hits=syn_dict[syn_dict.index.intersection(tok)]
		if len(hits)>0:
			g1g2= list(set(tok) - set(hits.index))
			#
			if len(abb.keys()) > 0:
				sus=set(set(hits.index) & abb.keys())
				#
				if len(sus) > 0:
					#
					def pw2(abb,al):
						s=pairwise2.align.globalms(abb, al,1, -1, -1, -0.5,score_only=True)
						s=s/len(abb)
						return s
					#
					# inspect abbreviation
					for a in sus:
						seq=abb[a]
						if type(syn_dict[a]) == str: # no paralogues
							sym=str(syn_dict[a])
							als=list(syn_dict[syn_dict == sym].index)
						else: # paralogues
							sym=syn_dict[a].tolist()
							als=set(itertools.chain.from_iterable([list(syn_dict[syn_dict == s].index) for s in sym]))
						#print(list(itertools.product([seq],als)))
						#try:
						if (m := max(starmap(pw2,itertools.product([seq],als)))) < 0.15:
							hits.drop(labels = [a],inplace=True)
							print("suspicious token/abbreviation:",a,abb[a],m)
						#else: 
							# remove sub-tokens of the abbreviation??? (may increase false negatives, but think of "IL2 RECEPTOR") 
						#
						#
						# except:
						# 	print("MAX ERRROR?")
						# 	print(idd)
						# 	print(a)
						# 	#
						# 	try:
						# 		print(list(starmap(pw2,itertools.product([seq],als))))
						# 	except:
						# 		print('LIST ERROR?')
						# 		print([seq])
						# 		print(als)
						# 		#print(pw2(*([seq],als])))
						# 	#
			hits=set(hits)
			g1g2=[w for w in g1g2 if re.match('[\w]+[\-\/][\w]+',w)]
			if len(g1g2) > 0:
				#	
				#print("A gene-gene occurrence(?)",g1g2)
				for pair in g1g2:
					pair=re.split('[\-\/]',pair)
					# BOTH must be gene aliases ("GENE-dirt" occurrences will be discarded)
					if all(p in syn_dict.index for p in pair):
						#
						hits.union(set(syn_dict[pair]))
			#
			hit=[(h,idd) for h in hits]
			return hit 
		#
		else:
			return []

	print("input %s corpus" % (tag))
	#
	#while (ids := next(chunks,False)):     
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
	maxl=len(list(pmid_to_abs.keys()))
	tot=0
	start=time.time()
	#
	# res=[]
	# inps=distro(list(pmid_to_abs.keys()),ncores,'abs_dict_','syn_dict_','abb_dict_')
	# while (inp := next(inps,False)):
	# 	print(inp)
	# 	res.append(list(itertools.starmap(find_in_series,inp)))
	# print(res)
	# res=[i for i in res if i]
	# print(res)
	# res=list(itertools.chain.from_iterable(res))
	# for t in res:
	# 	genes_pub[t[0]].append(t[1])

	# print("loop ends")
	# time.sleep(5)
	#
	with Pool(ncores) as pool:
		#res=[pool.apply(find_in_series,args=(idd,globals()["abs_dict_"+str(ids.index(idd))],globals()["syn_dict_"+str(ids.index(idd))])) for idd in ids]   
		print("starmap...")
		res=pool.starmap(find_in_series,distro(list(pmid_to_abs.keys()),ncores,'abs_dict_','syn_dict_','abb_dict_'))
		print("store in dictionary")
		res=[i for i in res if i]
		res=list(itertools.chain.from_iterable(res))
		for t in res:
			genes_pub[t[0]].append(t[1])
		#tot+=len(ids)/maxl
	#	
	#print("gene extraction", "%.2f percent done" % (100*tot))       
	#
	print('elapsed time, %d Abstracts : %f' % (maxl,time.time()-start))
	#
	genes=list(genes_pub.keys())
	#
	for key in genes:
		if genes_pub[key]==[]:
			del genes_pub[key]
	# update 
	genes=list(genes_pub.keys())
	#
	compress_pickle(wd+"data/"+tag+"_ids_to_genes",genes_pub)
	#
	## also save a list of all genes appearing at least once ## 
	tmp=re.sub(tag+"/$","",wd)
	with open(tmp+"NERed_Genes.txt","a") as nered:
		nered.write('\n'.join(genes))