# PMIDS to Abstract MESH
# exec(open("/home/musellla/tam_textmining/code/pmid_to_abs_mesh.py").read())
import os
import sys
#from nltk.corpus.reader.plaintext import PlaintextCorpusReader
import urllib
import urllib.request
import xmltodict
import dpath
#
def refparse(pmid):
	#
	url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&linkname=pubmed_pubmed_refs&id=%s" % (pmid)
	#
	file = urllib.request.urlopen(url)
	data = file.read()
	file.close()
	#
	data = xmltodict.parse(data)
	#print(data)
	refs=dpath.get(data,"**/Link") # list of {Id:value} dicts
	pmids=[ref['Id'] for ref in refs]
	#pmids=[int(p) for p in pmids]
	return pmids
#
refs=[]
print("retrieve references from %i articles..." % (len(sys.argv[2:])))
for p in sys.argv[2:]:
	refs+=refparse(p)
#
refs=list(set(refs))
#
print("%i PMIDs have been found" % (len(refs)))
#
with open(sys.argv[1]+".txt","w") as refs_file:
	refs_file.write('\n'.join(refs))
#
print("results stored in %s" % (sys.argv[1]+".txt"))
