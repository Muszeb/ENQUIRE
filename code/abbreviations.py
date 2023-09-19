abbreviations
### inspect abbreviations ###
sd=re.split("tmp",wd)[0]
#
aliasdf = pd.read_csv(sd+'input/aliases_uniprotensembl_as_df.tsv',sep='\t')
# remove 1-2 letters only aliases
aliases=aliasdf.alias_symbol.tolist()
aliases=[x for x in aliases if str(x) != 'nan']
#
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
#	
aliases=[a for a in aliases if a not in english_words_set]
aliases=[a for a in aliases if not bool(re.search("^[a-zA-Z]{1,2}$|^[\W\d]+$", a))]
aliases=[re.sub("[\(\)\[\]]",'',a) for a in aliases]
aliases=[a for a in aliases if a not in greek_alphabet.values()]
# remove aliases with greek letters not alphanumeric!
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
		mat=[i for i in mat if i !='']
		#
		if len(mat)>1:
			mat=[undash(w) for w in mat]
			mat=' '.join(mat)
			return ab+' '+mat+' '#+matout
		elif len(mat)==1:
			return mat[0]
		else:
			return None
#
def greek_al(v,find_greek=find_greek):
	v=v.split(';')
	v=[find_greek(a) for a in v]
	v=[a for a in v if a]
	return ';'.join(v)
aliases=[find_greek(a) for a in aliases]
aliases=[a for a in aliases if a]
#aliases=[re.sub(r'\W+', '', input) for input in aliases]
#aliases=[re.sub("[^\\w-\\s]+",'',a) for a in aliases]
# unless they are reference genes => regardless being a reference gene
#
def abbout(pmid,tokd=abs_dict,fuld=full_dict,abb=schwartz_hearst.extract_abbreviation_definition_pairs,aliases=aliases):
	#
	tok,ful=tokd[pmid],fuld[pmid]
	#
	abbs=abb(doc_text=ful)
	halt=list(set(tok) & set(abbs.keys())) # these are abbreviations to be checked
	# two-letters only are excluded anyways
	halt=[h for h in halt if len(h)>2]
	abbs=pd.Series(abbs)
	halt=abbs[halt].tolist()
	return halt
#
def greek_sub(m):
	if m.group() is not None:
		try: 
			return greek_alphabet[m.group().encode('utf-8')]
		except:
			return ''
#
halt=set(itertools.chain.from_iterable(map(abbout,list(abs_dict.keys()))))
halt=[re.sub('[\u0080-\uFFFF]',greek_sub,v) for v in halt]
aliases=[re.sub('[\u0080-\uFFFF]',greek_sub,v) for v in aliases]
	# apply pairwise2 to detect potential matches
def pw2(abb,al):
	s=pairwise2.align.globalms(abb, al, 1, -1, -.5, -.5,score_only=True)
	s=s/len(abb)
	return s
	
#

#
halt=halt+['Colony stimulating factor receptor 1', 'Colony stimulating factor', 'Hypoxia inducible protein','inducible nitric oxide synthase','tumor necrosis factor alpha']
halt=set(halt)
aliases=set(aliases)
print("compute scores for all halts")
matches=itertools.product(halt,aliases)
ncores=32
from multiprocessing import Pool
with Pool(ncores) as pool:
	#res=[pool.apply(find_in_series,args=(idd,globals()["abs_dict_"+str(ids.index(idd))],globals()["syn_dict_"+str(ids.index(idd))])) for idd in ids]   
	print("starmap...")
	res=list(pool.starmap(pw2,itertools.product(halt,aliases)))

#scores=list(itertools.starmap(pw2,matches))
#
import matplotlib.pyplot as plt # Set the figure size

plt.clf()
plt.rcParams["figure.figsize"] = [7.00, 3.50] 
plt.rcParams["figure.autolayout"] = True
# Data points for the histogram
bins=numpy.concatenate([numpy.arange(min(res),0,2.0),numpy.arange(0,max(res)+0.05,0.05)])
plt.hist(res,bins=bins) 
# Save the histogram 
plt.savefig('/home/musellla/histpos.png') # Display the plot plt.show()
#
prob=aliasdf[aliasdf['symbol']=='ARID3B']['alias_symbol'].tolist()
#
#AT-rich interactive domain-containing protein 3B
#Bdp
#
STOP 

# POSSIBILITY TO EXCLUDE ABBREVIATIONS

#pairs = schwartz_hearst.extract_abbreviation_definition_pairs(doc_text='The emergency room (ER) was busy')

################# code to use for converting pmid_to_kw format