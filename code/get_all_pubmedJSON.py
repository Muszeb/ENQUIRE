import requests
import json
import threading
from multiprocess import Process, Manager
import multiprocessing as mp
from multiprocessing import Pool
#
search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&mindate=1800/01/01&maxdate=2022/12/31&usehistory=y&retmode=json"
search_r = requests.post(search_url)
search_data = search_r.json()
webenv = search_data["esearchresult"]['webenv']
total_records = int(search_data["esearchresult"]['count'])
fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmax=9999&query_key=1&webenv="+webenv
#
ncores=4
chunksize=int(total_records/4)
chunks=list(range(0, total_records, chunksize))
#
def halfsearch(m):
	global chunks
	global_lock = threading.Lock()
	#print(len(nchunks))
	for i in range(chunks[m-1],chunks[m],10000):
		this_fetch = fetch_url+"&retstart="+str(i)
		print("Getting this URL: "+this_fetch)
		try: 
			fetch_r = requests.post(this_fetch)
			#print(fetch_r)
			if fetch_r.status_code==200:
				with global_lock:
					print("writing json/xml, chunk %s" % str(m))			
					with open("/home/musellla/pubmed_batch.json", 'a') as f:
						f.write(fetch_r.text)
			else:
				with global_lock:
					print("writing unretrieved batch id")	
					with open('/home/musellla/pubmed_unretrieved_batch.txt', 'a') as f:
						f.write(m+':'+i+'\n') # chunk:batch information
		except Exception as e:
			print("Error! Bad internet connection?")
			print(e)
			time.sleep(5)	
	   
with Pool(processes=ncores) as pool:
		#
		pool.map(halfsearch,list(range(1,len(chunks))))   
	
