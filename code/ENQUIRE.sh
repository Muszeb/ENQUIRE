#!/bin/bash
# txtmining_master (bash)
############################################################
# Help section                                             #
Help(){
	echo
	echo "####################################################################################"
	echo
	echo "Expanding Networks by Querying Unexpectedly Inter-Related Entities"
	echo
	echo "####################################################################################" 
	echo 
	echo "####################################################################################"
	echo
	echo "Usage: ./code/ENQUIRE.sh (-options...), " #" <working directory> <input PubMed_abstracts.txt>"
	echo
	echo -e "Legend:\t[-flag_short|--flag_long|config file variable, if available]:"
	echo
	#echo "<working directory> = "
	#echo "<input PubMed_abstracts.txt> =	"
	echo -e "[-p|--path|wd] = the path to the working directory (wd).\n\tIt must be the ENQUIRE main folder, with ./code and ./input as subfolders.\n\tThe default is the current working directory.\n"
	echo -e "[-i|--input|to_py] = input.txt: a 'seed' input text file containing one PMID per line.\n\tIt can be obtained from a PubMed querying specifying 'PMID' as the download format option.\n\tA minimun of 3 entries is required, but a list at least a few dozens articles is highly recommended.\n"
	echo -e "[-t|--tag|tag] = A characteristic tag definining the task.\n\tIt must be an alphanumeric string.\n"
	echo -e "[-j|--ncores|ncores] = The max number of CPU cores to be used.\n\tDefault is 6.\n"		
	echo -e "[-c|--combine-set|comb] = how many N entities to intersect to construct a query?\n\t3: loose searches, 4: moderate (default), 5: very strict queries.\n" 
	echo -e "[-r|--representativeness|thr] = representativeness threshold (%) for a subgraph to be included in the network expansion steps? (default: 0 %).\n\tExample: if a subgraph contains nodes exclusively mentioned in 10 papers out of a total of 100, that subgraph has a 10% representativeness.\n" 
	echo -e "[-a|--attempts|A] = how many query attempts (i.e. pairs of motifs or genes) should be run in order to connect any two subgraphs?\n\t1: conservative, 2: moderate (default), 3: greedy.\n"
	echo -e "[-k|--connectivity|K] = minimal community connectivity (K), which applies to any expansion-derived entities:\n\teach gene/MeSH term must be connected to at least K original communities to be incorporated in the expanded network - default: 2.\n"
	echo -e "[-f|--config] = if a config file is being used, specify its full path (e.g. input/textmining_config.txt).\n\tThis option overwrites any parameter set by a different option.\n"
	echo -e "[-w|--rscript|rscript] = path to the Rscript compiler of an installed R version (Should be the one installed using 'input/ENQUIRE.yml').\n\tThe default is 'which Rscript').\n"
	echo -e "[-e|--entity|etype] = which entity type (gene/MeSH) are you interested into? Omit or 'all' to textmine both entities.\n"
	echo -e "[-h|--help] = print this help message."
	echo
	echo "You might be seeing this Help because of an input error."
	echo
	echo "####################################################################################"
	exit 0
}
############################################################


############################################################
# Main program                                             #
############################################################
############################################################

# deal with options first
# This ensures we are not contaminated by variables from the environment.

while getopts ":hp:i:t:j:c:r:a:k:w:f:e::" options;do
	case ${options} in
		h|--help)
			# echo
			# echo "####################################################################################"
			# echo
			# echo "Iterative correlate-query-expand network construction by literature texmining"
			# echo
			# echo "####################################################################################" 
			# echo 
			# echo "####################################################################################"
			# echo
			# echo "Syntax: txtmining_master.sh (-options...), " #" <working directory> <input PubMed_abstracts.txt>"
			# echo
			# echo "Where"
			# echo
			# #echo "<working directory> = "
			# #echo "<input PubMed_abstracts.txt> =	"
			# echo "[-p|--path] = the path to the working directory (wd) : it's supposed to be the git repo, with ./code and ./input as subfolders"
			# echo "[-i|--input] = input.txt (to be saved under "wd/input/"): a text file containing one PMID per line. It can be obtained from PubMed advanced queries specifying "PMID" as the downloading option"
			# echo "[-t|--tag] = A characteristic tag definining the task"
			# echo "[-j|--ncores] = The max number of CPU cores to be used"
			# echo "[-s|--server] = which server are you running the script from? (options: srv1, srv3)"			
			# echo "[-c|--combine] = how many N entities to intersect to construct a query? (4: default and recommended, 3: loose searches, 5:very strict queries)" 
			# echo "[-r|--representativeness] = representativeness threshold (%) for a subgraph to be considered for querying? (default: 0 %). e.g. if 10 papers out of 100 contain nodes belonging to a subgraph, that subgraph has 10% representativeness" 
			# echo "[-a|--attempts] = how many query attempts (i.e. pairs of motifs or genes) should be run in order to connect any two subgraphs? (1: conservative, 2: moderate, 3: greedy)"
			# echo "[-k|--connectivity] = minimal connectivity (K), defined as the product of nodes and edges numbers, that all subgraphs are required to possess (applies to gene-only network expansion). Examples: K=2 for 2 nodes / 1 edge, K=6 for 3 nodes / 2 edges, K=9 for 3 nodes / 3 edges"
			# echo "[-f|--config] = if a config file is being used, specify its full path"
			# echo "[-w|--password] = password to access MariaDB SQL database of PubMed literature (set only if -f option have been chosed)"
			# echo "[-h|--help] = print help"
			# echo
			# echo "you might be seeing this Help because of an input error."
			# echo
			# echo "####################################################################################"
			# exit 1
			Help
			;;
		i | --input) 
			input=$OPTARG
			to_py=$input
			echo "input set to $input"
			;;
		p | --path)
			wd=$OPTARG
			cd "$wd"
			echo "working directory set to $wd"
			;;
		t | --tag) 
			tag=$OPTARG
			echo "tag set to $tag"
			;;
		j | --cores) 
			ncores=$OPTARG
			echo "$ncores CPU cores will be used"
			;;	
		c | --combine) 
			comb=$OPTARG
			echo "$comb entities out of a motifs pairs will generate an esearch query"
			;;
		r | --representativeness) 
			thr=$OPTARG
			echo "threshold for motif representativeness set to $thr percent"
			;;
		a | --attempts) 
			A=$OPTARG
			echo "$A query attempts for joining two subgraphs"
			;;
		k | --subgraph_connectivity) 
			K=$OPTARG
			echo "$K defines the subgraph connectivity"
			;;
		w | --rscript) 
			echo "set Rscript compiler"
			rscript=$OPTARG			
			;;	
		f | --config) 
			config=$OPTARG
			echo "config file with path $config has been passed"
			source "$config"
			;;
		e | --entity) 
			etype=$OPTARG
			if [[ -n $etype ]];then
				echo "Entity type will be restricted to $etype"
			fi
			;;							
		:)
			echo "Invalid option: $OPTARG requires an argument" 1>&2
			Help
			exit 1
			;;
		\?)
			echo "Invalid option: $OPTARG" 1>&2
			Help
			exit 1
			;;
	esac
done

### WARNING/ERROR HANDLING ###

if [[ -n "$config" ]]; then
	if [[ "$(grep -c "^[^\-\#]*=[^\-\#]*$\|^$" "$config")" -lt "$(grep -c ".*" "$config")" ]]; then
		echo "ERROR: potentially malicious config file detected!"
		echo "Only lines implying 'variable=value' are accepted!"
		echo -e "Your Config file:\n\n$(cat "$config")\n"
		Help
		exit 1
		### rscript required ###
	else
		echo -e "sourcing parameters from ${config}"
		source "$config"
		
	fi
fi

if [[ -z "$wd" ]]; then
	echo -e "Warning:\nMissing working directory, ENQUIRE will use the current working directory"
	wd=$(pwd)/
fi

if [[ -z "$rscript" ]]; then
	echo -e "Warning:\nRscript path was not specified under the 'rscript' variable,\n set to environment-available R interpreter"
	rscript=$(which Rscript)
	echo "$rscript"
fi

# halting errors

if [[ ! -f "$rscript" ]];then
	echo "ERROR! Rscript path is unvalid:"
	echo "$rscript"
	Help
	exit 1	
elif [[ -z "$to_py" ]]; then
	echo "Error: missing an input (required)"
	Help
	exit 1
elif [[ -z "$tag" ]]; then
	echo "Error: missing a tag/job name (required)"
	Help
	exit 1
### BLOCK ANY INPUT THAT IS SMALLER THAN 3 PAPERS ###
elif [[ "$(wc -l < "$to_py")" -lt 3  ]]; then
	echo "STOP: LESS THAN 3 ARTICLES HAVE BEEN PASSED"
	echo "A minimum of 3 PMIDs is required, but a list of a few dozens entries is recommended"
	exit 1
else
	echo "the script will run"
	echo
	tmp="${wd}tmp-${tag}"
	echo
	echo "create a temporary directory"
	# save tmp in its own variabl
	mkdir -p "$tmp"
	echo "create record of already visited pairs of motifs"
	if [[ -f "${tmp}/new_iteration.txt" ]];then
		touch "${tmp}/visited_motifs.tsv"
	fi
	echo
	echo "set Rscript compiler - don't read .Rprofile"
	rscript="${rscript} --no-environ"
	#rscript="/home/musellla/miniconda3/envs/R363/lib/R/bin/Rscript"
	echo "create a copy of source PMIDs"
	#touch ${wd}tmp/source_pmids.txt
	cat "$to_py" > "${tmp}/source_pmids.txt"
	echo
	echo "create an efetch inputs directory"
	rm -rf "${tmp}/efetch_inputs"
	mkdir -p "${tmp}/efetch_inputs"
	cp -f "$to_py" "${tmp}/efetch_inputs"
	to_py="$(basename $to_py)"
	#
	touch "${tmp}/previous_iteration.txt"
	# check if an expansion has already happened, hence update tag #
	it=0
	i=0
	#
	if [[ -f "${tmp}/new_iteration.txt" ]];then
		tag=$(cat "${tmp}/new_iteration.txt") # current
		prev=$(cat "${tmp}/previous_iteration.txt")
		# then also update it and i
		i=$(($(find "${tmp}/" -maxdepth 1 -type d | wc -l)-2)) # -2 because of "." and "efetch inputs"
		it=$i
	fi
	### DEFAULT VALUES ###
	if [[ -z "$comb" ]];then
		comb=4
		echo "default: ${comb} entities out of a motifs pairs will generate an esearch query"
	fi
	if [[ -z "$thr" ]];then
		thr=1
		echo "default: representativeness threshold set to $thr"	
	fi
	if [[ -z "$A" ]];then
		A=2
		echo "default: ${A} attempt per pair of subgraphs will be performed"	
	fi
	if [[ -z "$etype" ]];then
		etype='all'
		echo "default: both genes and MeSH terms will be collected"	
	fi
	if [[ -z "$rscript" ]];then
		rscript=$(which Rscript)
		echo "default: ${rscript} will be used as R interpreter"	
	fi
	if [[ -z "$ncores" ]];then
		ncores=6
		echo "default: ${ncores} cores will be used"	
	fi
	if [[ -z "$K" ]];then
		K=3
		echo "default: connectivity requirement set to ${K} communities"	
	fi
	#echo "number of iterations set to $iter"
	# status.log controls number of iterations
	rm -f -- "${tmp}/${tag}/status.log"
	mkdir -p "${tmp}/${tag}/data"
	mkdir -p "${tmp}/${tag}/data/subgraphs"
	touch "${tmp}/${tag}/status.log"
	# initiate query record #
	echo -e "PMID\tQuerySet" > "${tmp}/efetch_inputs/QueryToPMID_record.tsv"
	paste <(tr -d '\r' < "${tmp}/efetch_inputs/${to_py}") \
	 <(printf "${to_py}\n%.0s" $(seq 1 "$(nl "${tmp}/efetch_inputs/${to_py}" | wc -l)")) >> "${tmp}/efetch_inputs/QueryToPMID_record.tsv"
	### set maximum stack size (expressed in Kb)
	#ulimit -s 33554432
	### SAVE INPUT PARAMETERS! ###
	echo "Save Parameters in ad hoc file"
	echo -e "wd=${wd}\ntag=${tag}\nncores=${ncores}\nK=${K}\nA=${A}\nthr=${thr}\ncomb=${comb}\nto_py=${input}\netype=${etype}\nrscript=${rscript}" > "${tmp}/ENQUIRE_input_parameters.txt"
	echo
	###
	while  [[ "$(wc -w < "${tmp}/${tag}/status.log")" -lt 1 ]];
		do			
			echo
			echo "create a data directory and subgraphs subdirectory, if non existing"
			mkdir -p "${tmp}/${tag}/data"
			mkdir -p "${tmp}/${tag}/data/subgraphs"
			echo
			echo "##############################################################"
			echo
			echo "Fetch MESH in PMIDs and extract abstract"
			echo
			if [[ ! -f "${tmp}/${tag}/data/${tag}_ids_to_full.pbz2" ]]; then			
				#if python $wd"code/pmid_to_abs_mesh.py" ${tmp}/${tag}/ $tag "${tmp}/efetch_inputs/${to_py}"; then
		    	if python "${wd}/code/EDirect_pmid_to_abs_mesh.py" "${tmp}/${tag}/" "$tag" "${tmp}/efetch_inputs/${to_py}" "$ncores" "$etype" > "${tmp}/${tag}/input.log"; then
		 			echo "$(cat "${tmp}/${tag}/input.log")"
		 			rm -f "${tmp}/${tag}/input.log"   	
		    		echo "Exit code of 0, success"
				elif [[ "$(cat "${tmp}/${tag}/input.log")" = "ERROR: INVALID INPUT FILE" ]];then
					echo "FATAL: your input looks like:"
					echo "$(head -10 "${tmp}/efetch_inputs/${to_py}")"
					echo -e "\n\n"
					echo "while it should look like:"
					echo -e "\n2804834\n3580523\n3725990\n479723\n4528080\n"
					exit 1
				else
		    		echo "Exit code of $?, failure"
		    		set -e
		    		exit 1
				fi
			fi
			
			echo 
			echo "Completed - PMIDs_to_MESH and PMIDs_to_Abstract objects stored"
			echo
			echo "##############################################################"
			echo
			echo "##############################################################"
			echo
			echo "Tokenization: from full abstracts texts to lists of entities"
			echo
			if [[ ! -f "${tmp}/${tag}/data/${tag}_ids_to_abs.pbz2" ]]; then			
				if python "${wd}/code/abs_tokenize.py" "${tmp}/${tag}/" "$tag" "$ncores" "$etype"; then
		 			echo "Exit code of 0, success"
				else
		    		echo "Exit code of $?, failure"
		    		set -e
		    		exit 1
				fi
			fi
			echo 
			echo "Completed - Tokens generated"
			echo
			echo "##############################################################"
			echo
			echo "##############################################################"
			echo
			echo "Gene extraction from abstracts"
			echo

			if [[ ! -f "${tmp}/${tag}/data/${tag}_ids_to_genes.pbz2" ]]; then
				if python "${wd}/code/gene_extraction_v3.py" "${tmp}/${tag}/" "$tag" "$ncores"; then
		    		echo "Exit code of 0, success"
				else
		    		echo "Exit code of $?, failure"
		    		set -e
		    		exit 1
				fi
			fi
			echo 
			echo "Completed - PMIDs_to_Genes object stored"
			echo
			echo "##############################################################"
			echo
			echo "##############################################################"
			echo
			echo "Fisher Test for co-occurrence of entities in a contingency table"
			echo
			if [[ ! -f ${tmp}/${tag}/data/${tag}_allxall_sig_combinations_bh_to_cytoscape.tsv ]]; then
				echo "Remove any TSVs containing Fisher test statistics, if any"
				rm -f -- "${tmp}/${tag}/data/$tag_"*"allxall"*".tsv" 
				#mkdir -p "${wd}${tag}/data/subgraphs/"				
				if python "${wd}code/create_network_from_dict_URNV2.py" "${tmp}/${tag}/" "$tag" "$ncores" "$i"; then
		    		echo "Exit code of 0, success"
		    	elif (($? == 11)); then
		    		echo "STOP: no edges found in the graph" > "${tmp}/${tag}/status.log"
		    		break
				else
		    		echo "Exit code of $?, failure"
		    		set -e
		    		exit 1
				fi
			fi
			# CHECK IF THERE ARE ANY SIG CO-OCCURRENCES
			if [[ ! -f "${tmp}/${tag}/data/${tag}_allxall_sig_combinations_bh_to_cytoscape.tsv" ]] && [[ "$it" -lt 1 ]]; then
				echo "Error: No significantly co-occurring entities found" > "${tmp}/${tag}/status.log"
				break
			fi
			echo "Completed - Test statistics, adjustments edge and node tables created"
			echo
			echo "##############################################################"
			echo
			echo "##############################################################"
			echo
			# ALTERNATIVE to code at lines 323-334

			if [ "$it" != 0 ]; then
				echo "concatenate edges and nodes from previous iterations..."
				echo 
				if csvstack "${tmp}/${prev}/data/${prev}_allxall_sig_combinations_bh_to_cytoscape.tsv" "${tmp}/${tag}/data/${tag}_allxall_sig_combinations_bh_to_cytoscape.tsv" > "${tmp}/${tag}/data/tmp.tsv"; then
					mv -bf "${tmp}/${tag}/data/tmp.tsv" "${tmp}/${tag}/data/${tag}_allxall_sig_combinations_bh_to_cytoscape.tsv" 
					echo "edges ok"

				else
					echo "Exit code of $?, failure"
					set -e
					exit 1
				fi
				#
				if csvstack "${tmp}/${prev}/data/${prev}_allxall_pre_nodes_table_bh_to_cytoscape.tsv" "${tmp}/${tag}/data/${tag}_allxall_pre_nodes_table_bh_to_cytoscape.tsv" > "${tmp}/${tag}/data/tmp.tsv"; then
					mv -bf "${tmp}/${tag}/data/tmp.tsv" "${tmp}/${tag}/data/${tag}_allxall_pre_nodes_table_bh_to_cytoscape.tsv"
					echo "nodes ok"
				else
					echo "Exit code of $?, failure"
					set -e
					exit 1
				fi
			fi

			echo
			echo "Assign network-based parameters to nodes"
			if [[ ! -f "${tmp}/${tag}/data/${tag}_allxall_nodes_table_bh_to_cytoscape.tsv" ]]; then
				if python "${wd}code/fisher_to_netStatsV4.py" "${tmp}/${tag}/" "$tag" 3 "$thr"; then
		    		echo "Exit code of 0, success"
				else
		    		echo "Exit code of $?, failure"
		    		set -e
		    		exit 1
				fi
			fi	
			#
			echo "Network-based Nodes Prioritization (SANTA)"
			# check if subgraphs have been already ranked
			cols=$(head -n1 "${tmp}/${tag}/data/${tag}_allxall_nodes_table_bh_to_cytoscape.tsv" | wc -w)
			if [[  "$cols" -lt 7 ]]; then			
				if $rscript "${wd}/code/SANTA_NetRankings.R" "${tmp}/${tag}/" "$tag" incomplete "$A"; then
		    		echo "Exit code of 0, success"
				else
		    		echo "Exit code of $?, failure"
		    		set -e
		    		exit 1
				fi
				# if ${wd}/code/SANTA_subgraphs_${srv}.R ${tmp}/${tag}/ $tag $A; then
		  #   		echo "Exit code of 0, success"
				# else
		  #   		echo "Exit code of $?, failure"
		  #   		set -e
		  #   		exit 1
				# fi
			fi
			echo 
			echo "Completed - Rankings computed"
			echo
			echo "##############################################################"
			echo
			echo
			echo "##############################################################"
			echo
			echo "Gene-MeSH-Gene motives extraction"
			echo
			# x catches the event of a connected gene-mesh graph, thus triggering
			# the "merge gene-only graphs" section of the pipeline
			rm -f -- "${tmp}/${tag}/status.log"
			#
			python "${wd}/code/networkx_find_sub_motifs.py" "${tmp}/${tag}/" "$tag" "$thr" "$A" 2> "${tmp}/${tag}/status.log"
			echo "##############################################################"
			echo 
			it=$((it + 1))
			#rm -rf -- "${tmp}/${tag}/data/"*".pbz2"
			x="Gene-MeSH Network Completed - start connecting gene-only subgraphs"
			if ! [ "$(cat "${tmp}/${tag}/status.log")" = "$x" ]; then
				echo
				echo "this is status.log: $(cat "${tmp}/${tag}/status.log")"
				echo
				echo "##############################################################"
				echo

				echo "Reconstruct Cliques-based Network Communities"
				echo
				if ! ([ -f "${tmp}/${tag}/${tag}_Complete_edges_table_subgraph.tsv" ] || [ -f "${tmp}/${tag}/${tag}_ordered_queries.tsv" ]) \
					&& [ -f "${tmp}/${tag}/data/${tag}_allxall_sig_combinations_bh_to_cytoscape.tsv" ] ; then							
					echo "run incomplete version"
					if  $rscript "${wd}/code/SANTA_on_cliques-cliques_TSP.R" "${tmp}/${tag}/" "$tag" incomplete "$ncores" "$A" "$comb"; then
			    		echo "Exit code of 0, success"
			    	elif (($? == 12)); then
			    		echo "STOP: no cliques have been found" > "${tmp}/${tag}/status.log"
			    		break
			    	elif (($? == 13)); then
			    		echo "STOP: there is only one community" > "${tmp}/${tag}/status.log"
			    		break
					else
			    		echo "Exit code of $?, failure"
			    		set -e
			    		exit 1
					fi
				elif ! [ -f "${tmp}/${tag}/${tag}_ordered_queries.tsv" ] && \
					[ -f "${tmp}/${tag}/${tag}_Complete_edges_table_subgraph.tsv" ]; then							
					echo "run complete version"
					if $rscript "${wd}/code/SANTA_on_cliques-cliques_TSP.R" "${tmp}/${tag}/" "$tag" complete "$ncores" "$A" "$comb"; then
			    		echo "Exit code of 0, success"
			    	elif (($? == 12)); then
			    		echo " STOP: no cliques have been found" > "${tmp}/${tag}/status.log"
			    		break
			    	elif (($? == 13)); then
			    		echo "STOP: there is only one community" > "${tmp}/${tag}/status.log"
			    		break
			    	elif (($? == 14)); then
			    		echo "STOP: no significant clique clusters found" > "${tmp}/${tag}/status.log"
			    		break		
					else
			    		echo "Exit code of $?, failure"
			    		set -e
			    		exit 1
					fi
				else 
					echo "Found a pre-computed output"
					#set -e
					#exit 1
				fi
				echo
				echo "Completed - Candidate Queries to PubMed hav been generated"
				echo
				echo "#############################################################"
				echo
				echo "#############################################################"
				echo
				echo "Formulate PubMed Queries"
				echo
				if [[ ! "$tag" == "$(cat "${tmp}/previous_iteration.txt")" ]]; then		
					if python "${wd}code/interquery_from_clique_stats_TSP_EDirect.py" "${tmp}/${tag}/" "$tag" "${tmp}/source_pmids.txt" "$i" "${tmp}/${tag}/${tag}_ordered_queries.tsv" "$comb" "$A"; then
			    		echo "Exit code of 0, success"
					else
			    		echo "Exit code of $?, failure"
			    		set -e
			    		exit 1
					fi				
				fi
				echo
				echo -e "Completed - A New Set of $(grep -c . "${tag}_pmids_to_efetch.txt") PMIDs have been created"
				echo
				echo "##############################################################"	
				echo
				# keep track of new and old tags for merging dicts later on
				tag=$(cat "${tmp}/new_iteration.txt") # current
				prev=$(cat "${tmp}/previous_iteration.txt")
				to_py="${tag}_pmids_to_efetch.txt" 
				echo "Iteration number $it completed!"
				echo "next is $tag"
				echo "next PMID list is $to_py"
				echo
				i=$((i + 1))
				#it=$((it + 1))
				# check if file is not empty 
				#lines=$(wc -l "${tmp}/efetch_inputs/${tag}_pmids_to_efetch.txt" | cut -f 1 -d " ")
				#[[ "$lines" -gt 0 ]] && echo "Completed - Ranked list of queries defined" || (echo -e "######## \n STOP: \n ${wd}input/${tag}_pmids_to_efetch.txt is empty \n ########" && exit 1)
				echo "clean directory..."
				mkdir -p "${tmp}/${tag}" 
				touch "${tmp}/${tag}/status.log"
				rm -f -- "${tmp}/${prev}/status.log"	
				#
				if [[ "$(wc -l < "${tmp}/efetch_inputs/${to_py}")" -lt 1 ]]; then
					echo -e "##################\nSTOP:\n ${tmp}/efetch_inputs/${to_py} is empty\n##################"
					mkdir -p "${tmp}/${tag}"
					echo "The algorithm stops, as no new PMID has been found" > "${tmp}/${tag}/status.log"
					echo -e "It was not possible to continue the network expasion, but you can check for interesting motifs in \n ${tmp}/${prev}/data/${prev}_interactive_Cliques_Network.html"
    				#echo "Exit code of 0, success"
    				break
				fi
			else
				echo "CHECKPOINT (2): ARE THERE ANY DIRECT GENE-GENE CO-OCCURRENCES?"
				echo
				if [[ "$(wc -l < "${tmp}/${tag}/${tag}_Genes_edges_table_subgraph.tsv")" -lt 2 ]]; then 
					echo "STOP: THERE ARE NO DIRECT CO-OCCURRENCES AMONG GENES,"
					echo "THEREFORE A SUBGRAPHS EXPANSION CANNOT BE PERFORMED"
					echo "The algorithm stops, as no gene-MeSH-gene motif can be found" > "${tmp}/${tag}/status.log"
					echo 
					break
				else
					echo "CHECKPOINT (2): PASS" 
				fi
			fi
			echo "##############################################################"
		#	
		done
	
	#################
	
	#################
	
	#################
	
	if [[ "$(cat "${tmp}/${tag}/status.log")" = "$x" ]]; then
		echo
		i=$((i + 1)) # impose merging of previous co-occurrences
		echo
		echo "###########################################################"
		echo
		echo "$x"
		echo
		echo "MERGE GENE SUBGRAPHS"
		echo 
		mkdir -p "${tmp}/${tag}/data/gene-subgraphs"
		echo
		python "${wd}/code/networkx_sub_gene_graphs.py" "${tmp}/${tag}/" "$tag" "$thr" "$K"
		echo
		while [[ "$(find "${tmp}/${tag}/data/gene-subgraphs" -type f | wc -l)" -gt 2 ]]; # 2 because one edge table and one node table
			do 
				echo "Network-based Nodes Prioritization (SANTA)"
				cols=$(head -n1 "${tmp}/${tag}/data/gene-subgraphs/${tag}_nodes_table_gene-subgraph_0.tsv" | wc -w)
				if [[  "$cols" -lt 10 ]]; then							
					echo
					if $rscript "${wd}/code/SANTA_NetRankings.R" "${tmp}/${tag}/" "$tag" "$A"; then
						echo
						echo "Exit code of 0, success"
					else
			    		echo "Exit code of $?, failure"
			    		set -e
			    		exit 1
					fi
				fi
				echo "Completed - Rankings computed"
				echo
				echo "############################################################"
				echo
				echo "############################################################"
				echo
				echo "Construct subgraphs-connecting esearch queries"
				echo "tag is $tag"
				next="${tag}_subgraphs_expansion"
				echo
				echo "##############################################################"
				echo
				echo "Reconstruct Cliques-based Network Communities"
				echo
				if ! ([ -f "${tmp}/${tag}/${tag}_Complete_edges_table_subgraph.tsv" ] || [ -f "${tmp}/${tag}/${tag}_ordered_queries.tsv" ]) \
					&& [ -f "${tmp}/${tag}/data/${tag}_allxall_sig_combinations_bh_to_cytoscape.tsv" ] ; then							
					echo "run incomplete version"
					if $rscript "${wd}code/SANTA_on_cliques-cliques_TSP.R" "${tmp}/${tag}/" "$tag" incomplete "$ncores" "$A" "$comb"; then
			    		echo "Exit code of 0, success"
			    	elif (($? == 12)); then
			    		echo "STOP: no cliques have been found" > "${tmp}/${tag}/status.log"
			    		break
			    	elif (($? == 13)); then
			    		echo "STOP: there is only one community" > "${tmp}/${tag}/status.log"
			    		break
					else
			    		echo "Exit code of $?, failure"
			    		set -e
			    		exit 1
					fi
				elif ! [ -f "${tmp}/${tag}/${tag}_ordered_queries.tsv" ] && \
					[ -f "${tmp}/${tag}/${tag}_Complete_edges_table_subgraph.tsv" ]; then							
					echo "run complete version"
					if $rscript "${wd}code/SANTA_on_cliques-cliques_TSP.R" "${tmp}/${tag}/" "$tag" complete "$ncores" "$A" "$comb"; then
			    		echo "Exit code of 0, success"
			    	elif (($? == 12)); then
			    		echo " STOP: no cliques have been found" > "${tmp}/${tag}/status.log"
			    		break
			    	elif (($? == 13)); then
			    		echo "STOP: there is only one community" > "${tmp}/${tag}/status.log"
			    		break
			    	elif (($? == 14)); then
			    		echo "STOP: no significant clique clusters found" > "${tmp}/${tag}/status.log"
			    		break				    		
					else
			    		echo "Exit code of $?, failure"
			    		set -e
			    		exit 1
					fi
				else 
					echo "Found a pre-computed output"
					#set -e
					#exit 1
				fi
				echo
				echo "Completed - Candidate Queries to PubMed hav been generated"
				echo
				echo "#############################################################"
				echo
				echo "#############################################################"
				echo
				echo "Formulate PubMed Queries"
				echo
				if [[ ! "$tag" == "$(cat "${tmp}/previous_iteration.txt")" ]]; then		
					if python "${wd}/code/interquery_from_clique_stats_TSP_EDirect.py" "${tmp}/${tag}/" "$tag" "${tmp}/source_pmids.txt" "$i" "${tmp}/${tag}/${tag}_ordered_queries.tsv" "$comb" "$A" ; then
			    		echo "Exit code of 0, success"
					else
			    		echo "Exit code of $?, failure"
			    		set -e
			    		exit 1
					fi				
				fi
				echo
				echo -e "Completed - A New Set of $(grep -c . "${tag}_pmids_to_efetch.txt") PMIDs have been created"
				echo
				echo "##############################################################"	
				# keep track of new and old tags for merging dicts later on
				prev=$(cat "${tmp}/previous_iteration.txt") 
				tag=$(cat "${tmp}/new_iteration.txt") # current
				#$(cat "${tmp}/previous_iteration.txt")
				to_py="${tag}_pmids_to_efetch.txt" 
				echo
				#echo "Iteration number $it completed!"
				echo "next is $tag"
				echo "next PMID list is $to_py"
				echo
				if [[ "$(wc -l < "${tmp}/efetch_inputs/${to_py}")" -lt 1 ]]; then
					echo -e "##################\nSTOP:\n${tmp}/efetch_inputs/${to_py} is empty\n##################"
					mkdir -p "${tmp}/${tag}"
					echo "The algorithm stops, as no new PMID has been found" > "${tmp}/${tag}/status.log"
					echo -e "It was not possible to continue the network expasion, but you can check for interesting cliques in \n ${tmp}/${prev}/${prev}_interactive_Cliques_Network.html"
    				#echo "Exit code of 0, success"
    				break
				fi
				echo
				echo "###########################################################"
				echo
				echo "##############################################################"
				echo
				echo
				echo "create a data directory and subgraphs subdirectory, if non existing"
				mkdir -p "${tmp}/${tag}/data"
				mkdir -p "${tmp}/${tag}/data/subgraphs"
				mkdir -p "${tmp}/${tag}/data/gene-subgraphs"
				echo
				#i=$((i + 1))
				it=$((it + 1))
				echo "##############################################################"
				echo
				echo "Fetch MESH in PMIDs and extract abstract"
				echo
				if [[ ! -f "${tmp}/${tag}/data/${tag}_ids_to_abs.pbz2" ]]; then			
					if python "${wd}/code/EDirect_pmid_to_abs_mesh.py" "${tmp}/${tag}/" "$tag" "${tmp}/efetch_inputs/${to_py}" "$ncores" "$etype"; then
			    		echo "Exit code of 0, success"
					else
			    		echo "Exit code of $?, failure"
			    		set -e
					fi
				fi			
				echo 
				echo "Completed - PMIDs_to_MESH and PMIDs_to_Abstract objects stored"
				echo
				echo "##############################################################"
				echo
				echo "##############################################################"
				echo
				echo "Tokenization: from full abstracts texts to lists of entities"
				echo
				if [[ ! -f "${tmp}/${tag}/data/${tag}_ids_to_abs.pbz2" ]]; then			
					if python "${wd}/code/abs_tokenize.py" "${tmp}/${tag}/" "$tag" "$ncores"; then
			 			echo "Exit code of 0, success"
					else
			    		echo "Exit code of $?, failure"
			    		set -e
			    		exit 1
					fi
				fi
				echo 
				echo "Completed - Tokens generated"
				echo
				echo "##############################################################"
				echo
				echo "##############################################################"
				echo
				echo "Gene extraction from abstracts"
				echo
				if [[ ! -f "${tmp}/${tag}/data/${tag}_ids_to_genes.pbz2" ]]; then
					if python "${wd}/code/gene_extraction_v3.py" "${tmp}/${tag}/" "$tag" "$ncores"; then
			    		echo "Exit code of 0, success"
					else
			    		echo "Exit code of $?, failure"
			    		set -e
			    		exit 1
					fi
				fi
				echo 
				echo "Completed - PMIDs_to_Genes object stored"
				echo
				echo "##############################################################"
				echo
				echo "##############################################################"
				echo
				echo "Fisher Test for co-occurrence of entities in a contingency table"
				echo
				if [[ ! -f "${tmp}/${tag}/data/${tag}_allxall_sig_combinations_bh_to_cytoscape.tsv" ]]; then
					echo "Remove any TSVs containing Fisher test statistics, if any"
					rm -f -- "${tmp}/${tag}/data/$tag_"*"allxall"*".tsv" 
					#mkdir -p "${wd}${tag}/data/subgraphs/"
					
					if python "${wd}code/create_network_from_dict_URNV2.py" "${tmp}/${tag}/" "$tag" "$ncores" "$i"; then
			    		echo "Exit code of 0, success"
			    	elif (($? == 11)); then
		    			echo "STOP: no edges found in the graph" > "${tmp}/${tag}/status.log"
		    			break				
					else
			    		echo "Exit code of $?, failure"
			    		set -e
			    		exit 1
					fi
				fi
				#
				if [ "$it" != 0 ]; then
					echo "concatenate edges and nodes from previous iterations..."
					echo 
					if csvstack "${tmp}/${prev}/data/${prev}_allxall_sig_combinations_bh_to_cytoscape.tsv" "${tmp}/${tag}/data/${tag}_allxall_sig_combinations_bh_to_cytoscape.tsv" > "${tmp}/${tag}/data/tmp.tsv"; then
						mv -bf "${tmp}/${tag}/data/tmp.tsv" "${tmp}/${tag}/data/${tag}_allxall_sig_combinations_bh_to_cytoscape.tsv" 
						echo "edges ok"

					else
						echo "Exit code of $?, failure"
						set -e
						exit 1
					fi
					#
					if csvstack "${tmp}/${prev}/data/${prev}_allxall_pre_nodes_table_bh_to_cytoscape.tsv" "${tmp}/${tag}/data/${tag}_allxall_pre_nodes_table_bh_to_cytoscape.tsv" > "${tmp}/${tag}/data/tmp.tsv"; then
						mv -bf "${tmp}/${tag}/data/tmp.tsv" "${tmp}/${tag}/data/${tag}_allxall_pre_nodes_table_bh_to_cytoscape.tsv"
						echo "nodes ok"
					else
						echo "Exit code of $?, failure"
						set -e
						exit 1
					fi
				fi
				#
				echo "Assign network-based parameters to nodes"
				if [[ ! -f "${tmp}/${tag}/data/${tag}_allxall_nodes_table_bh_to_cytoscape.tsv" ]]; then
					if python "${wd}code/fisher_to_netStatsV4.py" "${tmp}/${tag}/" "$tag" 3 "$thr"; then
		    			echo "Exit code of 0, success"
					else
		    			echo "Exit code of $?, failure"
		    			set -e
		    			exit 1
					fi
				fi	
				echo "Completed - Test statistics, adjustments edge and node tables created"
				echo
				echo "##############################################################"
				echo
				echo "##############################################################"
				echo
				echo "Tie core and expansion gene sets by shortest connecting paths"
				echo
				if [[ ! -f "${tmp}/${tag}/${tag}_Genes_table_subgraph.tsv" ]]; then
					echo
					if python "${wd}code/clique_expansion.py" "${tmp}/${tag}/" "$tag" "$prev" "$K"; then
			    		echo "Exit code of 0, success"
					else
			    		echo "Exit code of $?, failure"
			    		set -e
			    		exit 1
					fi
				fi

				echo "Completed - expanded network created (iteration ${i})"
				echo
				echo "Check if any relevant subgraph remains..."
				echo
				if [ -z "$(ls -A "${tmp}/${tag}/data/gene-subgraphs")" ]; then
					echo
					if python "${wd}code/networkx_sub_gene_graphs.py" "${tmp}/${tag}/" "$tag" "$thr" "$K"; then
			    		echo "Exit code of 0, success"
					else
			    		echo "Exit code of $?, failure"
			    		set -e
			    		exit 1
					fi
				fi
				echo
				i=$((i + 1))
				echo "##############################################################"
			echo
			done
			echo "##############################################################"
			echo
			echo "GENERATE FINAL PLOTS"
			echo
			if [[ "$(\ls "${tmp}/${tag}/")" = status.log ]];then
				tag="$prev"
			fi
			if ! ([ -f "${tmp}/${tag}/${tag}_Complete_edges_table_subgraph.tsv" ] || [ -f "${tmp}/${tag}/${tag}_ordered_queries.tsv" ]) \
				&& [ -f "${tmp}/${tag}/data/${tag}_allxall_sig_combinations_bh_to_cytoscape.tsv" ] ; then							
				echo "run incomplete version"
				if $rscript "${wd}code/SANTA_fina_lplots.R" "${tmp}/${tag}/" "$tag" incomplete "$ncores" "$A" "$comb"; then
		    		echo "Exit code of 0, success"
				else
		    		echo "Exit code of $?, failure"
		    		set -e
		    		exit 1
				fi
			elif ! [ -f "${tmp}/${tag}/${tag}_ordered_queries.tsv" ] && \
				[ -f "${tmp}/${tag}/${tag}_Complete_edges_table_subgraph.tsv" ]; then							
				echo "run complete version"
				if $rscript "${wd}code/SANTA_final_plots.R" "${tmp}/${tag}/" "$tag" complete "$ncores" "$A" "$comb"; then
		    		echo "Exit code of 0, success"
				else
		    		echo "Exit code of $?, failure"
		    		set -e
		    		exit 1
				fi
			else 
				echo "Found a pre-computed output"
				#set -e
				#exit 1
			fi
			echo			

		###
		# EXCEPTIONS HANDLING
			
		# elif [[ "$(wc -w ${tmp}/${tag}/status.log | cut -f 1 -d " ")" = 0 ]]; then
  #   		echo "Exit code of 0, success"
  #   		# check line number of the motives list
		# 	lines=$(wc -l "${tmp}/${tag}/data/${tag}_motifs_subs_and_ranks.tsv" | cut -f 1 -d " ")
		# 	[[ "$lines" -gt 1 ]] && echo "Completed - Motifs identified" || (echo -e "######## \n STOP: \n ${wd}data/${tag}_motifs_strengthonly.tsv is empty \n ########" && exit 1)
		# 	echo
		elif [[ "$(cat "${tmp}/${tag}/status.log")" = "The algorithm stops, as no new PMID has been found" ]]; then
    		echo
    		echo -e "It was not possible to continue the network expasion, but you can check for interesting cliques in \n ${tmp}/${prev}/${prev}_interactive_Cliques_Network.html"
    		echo
    		echo "##############################################################"
			echo
			# echo "GENERATE FINAL PLOTS"
			# echo
			# if [[ "$(\ls ${tmp}/${tag}/)" = status.log ]];then
			# 	rm -rf -- ${tmp}/${tag}/
   #  			tag=$prev				
			# fi
			# if ! ([ -f ${tmp}/${tag}/${tag}_Complete_edges_table* ] || [ -f ${tmp}/${tag}/${tag}_ordered_queries.tsv ]) \
			# 	&& [ -f ${tmp}/${tag}/data/${tag}_allxall_sig_combinations_bh_to_cytoscape.tsv ] ; then							
			# 	echo "run incomplete version"
			# 	if $wd"code/SANTA_final_plots.R" ${tmp}/${tag}/ $tag incomplete $ncores $A $comb; then
		 #    		echo "Exit code of 0, success"
			# 	else
		 #    		echo "Exit code of $?, failure"
		 #    		set -e
		 #    		exit 1
			# 	fi
			# elif ! [ -f ${tmp}/${tag}/${tag}_ordered_queries.tsv ] && \
			# 	[ -f ${tmp}/${tag}/${tag}_Complete_edges_table* ]; then							
			# 	echo "run complete version"
			# 	if $wd"code/SANTA_final_plots.R" ${tmp}/${tag}/ $tag complete $ncores $A $comb; then
		 #    		echo "Exit code of 0, success"
			# 	else
		 #    		echo "Exit code of $?, failure"
		 #    		set -e
		 #    		exit 1
			# 	fi
			# else 
			# 	echo "Found a pre-computed output"
			# 	#set -e
			# 	#exit 1
			# fi
			# echo			
			# if ! [ -f ${tmp}/${tag}/${tag}_Complete_edges_table* ]; then							
			# 	echo -e "WARNING: a complete network could not be achieved,\n hence the final graphs are poorly connected,\n and are therefore marked with a '_disc_' suffix"
			# 	echo "tag is: $tag"
			# 	if python $wd"code/sort_network.py" ${tmp}/${tag}/ $tag 'incomplete'; then
		 #    		echo "Exit code of 0, success"
			# 	else
		 #    		echo "Exit code of $?, failure"
		 #    		set -e
		 #    		exit 1
			# 	fi
			# 	echo
			# fi
    		# check line number of thes motives list
			#lines=$(wc -l "${tmp}/${tag}/data/${tag}_motifs_subs_and_ranks.tsv" | cut -f 1 -d " ")
			#[[ "$lines" -gt 1 ]] && echo "Completed - Motifs identified" || (echo -e "######## \n STOP: \n ${wd}data/${tag}_motifs_strengthonly.tsv is empty \n ########" && exit 1)
			echo				
		elif [[ "$(cat "${tmp}/${tag}/status.log")" = "Error: Significant co-occurrences only yielded a sparse network" ]]; then
			echo
			cat "${tmp}/${tag}/status.log"
			echo
			echo -e "It was not possible to continue the network expasion, but you can check for significantly co-occurring terms in \n\n ${tmp}/${tag}/data/${tag}_allxall_nodes_table_bh_to_cytoscape.tsv \n\n and \n\n ${tmp}/${tag}/data/${tag}_allxall_sig_combinations_bh_to_cytoscape.tsv"
			echo
		elif [[ "$(cat "${tmp}/${tag}/status.log")" = "The algorithm stops, as no gene-MeSH-gene motif can be found" ]]; then
			echo
			cat "${tmp}/${tag}/status.log"
			echo
			echo -e "It was not possible to continue the network expasion, but you can check for significantly co-occurring terms in \n\n ${tmp}/${tag}/${tag}_Complete_edges_table_subgraphs.tsv \n\n and \n\n ${tmp}/${tag}/${tag}_Complete_nodes_table_subgraphs.tsv"
			echo
			echo "Additionally, Literature support with hyperlinks can be found in ${tmp}/${tag}/${tag}_Complete_edges_literature_links.tsv"
			echo
		elif [[ "$(cat "${tmp}/${tag}/status.log")" = "Error: No significantly co-occurring entities found" ]]; then
			echo
			cat "${tmp}/${tag}/status.log"
			echo
			echo -e "It was not possible to continue the network expasion.\nThe problem might have occurred because the set of PMIDs (unexpectedly?) consists of a random sample of genes and concepts.\nA further explaination might be that relatively few papers had both an abstract and a valid MESH list, which is required for the test statistics.\nIn case of complaints, please report your experience."
			echo
		elif [[ "$(cat "${tmp}/${tag}/status.log")" = "STOP: no edges found in the graph" ]]; then
			echo
			cat "${tmp}/${tag}/status.log"
			echo
			echo -e "There were no significant co-occurrences and it was not possible to create a network.\nYou can see the test statistics in \n\n ${tmp}/${tag}/data/${tag}_edgelist_allxall.tsv\n\nThe problem might have occurred because the set of PMIDs (unexpectedly?) consists of a random sample of genes and concepts.\nA further explaination might be that relatively few papers had both an abstract and a valid MESH list, which is required for the test statistics.\nIn case of complaints, please report your experience."
			echo
		elif [[ "$(cat "${tmp}/${tag}/status.log")" = "STOP: no cliques have been found" ]]; then
			echo
			cat "${tmp}/${tag}/status.log"
			echo
			echo -e "The algorithm stops.\nThe problem might have occurred because the set of PMIDs (unexpectedly?) is small or consists of a random sample of genes and concepts.\nA further explaination might be that relatively few papers had both an abstract and a valid MESH list, which is required for the test statistics.\nIn case of complaints, please report your experience."
			echo
		elif [[ "$(cat "${tmp}/${tag}/status.log")" = "STOP: there is only one community" ]]; then
			echo
			cat "${tmp}/${tag}/status.log"
			echo
			echo -e "The network can not be expanded any further as it is either densely connected or very small."
			echo
		elif [[ "$(cat "${tmp}/${tag}/status.log")" = "STOP: no significant clique clusters found" ]]; then
			echo
			cat "${tmp}/${tag}/status.log"
			echo
			echo -e "The algorithm stops.\nThe problem might have occurred because the set of PMIDs (unexpectedly?) consists of an overly connected sample of genes and concepts.\nA further explaination might be that relatively few papers had both an abstract and a valid MESH list, which is required for the test statistics.\nIn case of complaints, please report your experience."
			echo
		else
    		echo 
    		echo "unknown status file entry:" #¯\_(ツ)_/¯
    		cat "${tmp}/${tag}/status.log"
    		set -e
    		exit 1
		fi	
	echo "##############################################################"
	echo
	echo "clean directories..."
	cd "$tmp"
	rm -f -- *_iteration.txt
	find . -name '*.pbz2' -type f -delete
	for dir in $(find . -maxdepth 1 -mindepth 1 -type d);do
		if [[ "$(\ls "${dir:?}/")" = status.log ]];then
			rm -rf -- "${dir:?}/"
		fi
	done
	echo
	echo "edit headers..."
	new_header_edge="KW1\tKW2\tP.value\tadj.P.value\tFisherP.value\tztPois.cdf\tRaw.CO.Occurrence\tExpansion"
	sed -i -e "1s/.*/$(echo -e "$new_header_edge")/" */*_edges_table_subgraph.tsv
	new_header_node="KW\tCategory\tRepresentativeness\tRaw.Occurrence\tBetweenness.centrality\tCloseness.centrality\tStrength\tWscore\tKnode\tKnode.rank\tCommunityMembership"
	sed -i -e  "1s/.*/$(echo -e "$new_header_node")/" */*Complete_nodes_table_subgraph.tsv
	new_header_node=$(echo "$new_header_node" | sed 's/\\tCommunityMembership//')
	sed -i -e "1s/.*/$(echo -e "$new_header_node")/" */*Genes_nodes_table_subgraph.tsv	
	sed -i -e "1s/.*/$(echo -e "$new_header_node")/" */*MeSH_nodes_table_subgraph.tsv	
	new_header_links="PMID\tKW1\tKW2\tlink\tlink_explicit\tQuerySet"
	sed -i -e "1s/.*/$(echo -e $new_header_links)/" */*_Complete_edges_literature_links.tsv
	echo
	echo
	echo "JOB FINISHED!" 
	echo	
	echo 
	echo "Completed - Textmined network of genes and MESH terms created"
	echo
	echo "##############################################################"		
fi