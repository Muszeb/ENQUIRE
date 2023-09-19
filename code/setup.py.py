#exec(open("/home/musellla/tam_textmining/code/setup.py.py").read())
import os
import sys
from nltk.corpus.reader.plaintext import PlaintextCorpusReader
import numpy as np 
import pandas as pd
import itertools
from itertools import islice
import multiprocessing as mp
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
import csv
try:
    import _pickle as pickle
except ModuleNotFoundError:
    	import pickle
#
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

#
english_words_set= english_words_set.union(set(words.words()))
#
wd='/home/musellla/tam_textmining/tmp-vcells/vcells/'
sd= re.split("tmp",wd)[0]
tag='vcells'
tmp=re.split('\/'+tag+'\/',wd)[0]+'/'
#
Entrez.email='Luca.musella@uk-erlangen.de'
Entrez.api_key='b89451c47164b8705c76108eec6386658607'
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
aliasdf = pd.read_csv(sd+'input/aliases_uniprotensembl.tsv',sep='\t')
# remove 1-2 letters only aliases
aliases=aliasdf.aliases.tolist()
# unless they are reference genes => regardless being a reference gene
#refs=aliasdf.symbol.tolist()
def alclean(al):
	al=al.split(';')
	al=[a for a in al if a not in english_words_set]
	al=[a for a in al if not bool(re.search("^[a-zA-Z]{1,2}$|^\d+$", a))]
	al=[re.sub("[\(\)\[\]]",'',a) for a in al]	
	return ';'.join(al)
#
aliasdf['aliases']=aliasdf.aliases.apply(alclean,1)
#
#dictionarize
#
#bad_hits=list(set(filter(lambda x: re.match("^[a-zA-Z]{1,2}$",str(x)),aliases)))
#aliasdf=aliasdf[~aliasdf['aliases'].isin(bad_hits)]
#
syms=aliasdf.symbol.tolist()
### load data ###
pmid_to_abs = decompress_pickle(wd+'/data/'+tag+'_ids_to_abs.pbz2')
###
with gzip.open(tmp+'ids_to_abs_raw.json.gz', 'rt', encoding='UTF-8') as zipfile:
    pmids_to_raw_abs = json.load(zipfile)
#
# ACHTUNG! KEYS ARE NOW STRINGS
pmids_to_raw_abs={int(k):v for k,v in pmids_to_raw_abs.items()} 
###
print("PMIDs and alias DF loaded")
#
print("set multithreading")
##
print("duplicate PMIDs-abstracts and aliasases dictionary for parallelization...")

ncores= 4#int(sys.argv[3]) #32

var_names_abs=["abs_dict_"+str(i) for i in range(ncores)]
var_names_syn=["syn_dict_"+str(i) for i in range(ncores)]

for name in var_names_abs:
	globals()[name] = pmids_to_raw_abs	

for name in var_names_syn:
	globals()[name] = aliasd=pd.Series(aliasdf['aliases'].tolist(),index=aliasdf.symbol).to_dict()


#
# to avoid overtaking the max number of threads, split the list (feed as first "climbME" argument)
genes_in_abs=dict({sym:[] for sym in syms})
chunks = iter([list(pmids_to_raw_abs.keys())[x:x+ncores] for x in range(0, len(list(pmid_to_abs.keys())), ncores)])
idlist=[32208142,
32663198,
31920150,
32929201,
33229588,
32785653,
33469052,
32103760,
32313727,
32802195]
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

def score(txt,kw,mode):
	if mode=='fast': # score only, check ig
		return pairwise2.align.localms(txt,kw, 1, -2,-5,-5, score_only=True)
				# normalization: best score is 1 (exact match)
def greek_sub(m):
	if m.group() is not None:
		try: 
			return greek_alphabet[m.group().encode('utf-8')]
		except:
			return None

id_abs={k:re.sub('[\u0080-\uFFFF]',greek_sub,v) for k,v in id_abs.items()}
txt=id_abs[32929201]