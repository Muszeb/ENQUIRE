#!/bin/bash

echo -e "############# TURN ENQUIRE NETWORKS INTO KNOWLEDGE GRAPHS USING NEO4J - UTILITY SCRIPT ##############"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export SCRIPT_DIR
echo "Path to code: $SCRIPT_DIR"

############################################################
# Main program                                             #
############################################################
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
	echo "Usage: ENQUIRE2KG.sh [script_arguments]" #" <working directory> <input PubMed_abstracts.txt>"
	echo
	echo -e "Legend:\t[-flag_short|--flag_long|config file variable, if available]:"
	echo
	#echo "<working directory> = "
	#echo "<input PubMed_abstracts.txt> =	"
	echo -e "[-i|--image|image] = the path to the singularity image file (.sif). Defaults to 'ENQUIRE.sif'.\n"
	echo -e "[-p|--path|wd] = the path to the working directory (wd), where the output directory will be written in.\n\tIt must be the ENQUIRE main folder, with ./code and ./input as subfolders.\n\tThe default is the current working directory.\n"
	echo -e "[-t|--tag|tag] = A tag definining the task.\n\tIt must be an alphanumeric string (underline_spaced_words are accepted).\n"
	echo -e "[-d|--inputdir|input] = path to the input data folder. It must point to an ENQUIRE-generated directory containing co-occurrence network data\n\t(e.g https://github.com/Muszeb/ENQUIRE/tree/main/tmp-Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System).\n"
	echo -e "[-f|--config] = if a config file is being used, specify its full path (e.g. input/textmining_config.txt).\n\tThis option overwrites any parameter set by a different option.\n"
	echo -e "[-h|--help] = print this help message."
	echo
	echo "You might be seeing this Help because of an input error."
	echo
	echo "####################################################################################"
	exit 0
}

while getopts ":hi:p:t:f:d::" options;do
	case ${options} in
		h|--help)
			Help
			;;
		i | --image)
			sif=$OPTARG
			echo "SIF set to $wd"
			;;
		p | --path)
			wd=$OPTARG
			cd "$wd"
			echo "working directory set to $wd"
			;;
		d | --inputdir) 
			input=$OPTARG
			echo "input directory set to $input"
			;;
		t | --tag) 
			tag=$OPTARG
			echo "tag set to $tag"
			;;		
		f | --config) 
			config=$OPTARG
			echo "config file with path $config has been passed"
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

	if [[ -f "$config" ]]; then
		if grep -q $'\r' "$config" && ! grep -q $'\r\n' "$config"; then
			echo "File '$config' uses CR (Older Macintosh) line endings. Converting to LF (Unix)."
			sed -i 's/\r/\n/g' "$config"
		elif grep -q $'\r\n' "$config"; then
			echo "File '$config' uses CR LF (Windows) line endings. Converting to LF (Unix)."
			sed -i 's/\r$//' "$config"
		elif grep -q $'\n' "$config"; then
			echo "File '$config' uses LF (Unix) line endings."
		else
			echo "Unable to determine EOL style or file '$config' is empty."
		fi
	else
		echo "Config file '$config' not found or specified."
		exit 1
	fi

	if [[ "$(grep -c "^[^\-\#]*=[^\-\#]*$\|^$" "$config")" -lt "$(grep -c ".*" "$config")" ]]; then
		echo "ERROR: potentially malicious config file detected!"
		echo "Only lines implying 'variable=value' are accepted!"
		echo -e "Your Config file:\n\n$(cat "$config")\n"
		Help
		exit 1
	fi
	### rscript required ###
	echo -e "sourcing parameters from ${config}"
	source "$config"
fi

if [[ -z "$wd" ]]; then
	echo -e "Warning:\nMissing working directory, ENQUIRE will use the current working directory"
	wd=$(pwd)/
fi

if [[ -z "$sif" ]]; then
	echo -e "Warning:\nMissing SIF specification, will be set to 'ENQUIRE.sif' and assumed to be in the current working directory"
	sif="ENQUIRE.sif"
fi

# if [[ -z "$sd" ]]; then
# 	echo -e "Warning:\nMissing input data directory, ENQUIRE will use the current working directory"
# 	sd=$(pwd)/
# 	export sd
# fi

# halting errors
	
if [[ -z "$input" ]]; then
	echo "Error: missing an input (required)"
	Help
	exit 1
fi

if [[ -z "$tag" ]]; then
	echo "Error: missing a tag/job name (required)"
	Help
	exit 1
fi

###### ADDITIONAL PARAMETER SANIFICATION #######
if [[ "$tag"  =~ [^a-zA-Z0-9_] ]];then
	echo "FATAL: tag must be an alphanumeric string (underline-spaced words are accepted)"
	Help
	exit 1
fi

###### Passed halting errors ##### 
echo "the script will run"
echo
export tmp="${wd}enquire2kg-${tag}"
echo
echo "create a temporarily mounted directory"
# save tmp in its own variable

if [ -d "$tmp" ]; then
  echo "Warning: $tmp already exists. If containing a previously computed output, it may cause errors during knowledge graph generation."
fi

### SAVE INPUT PARAMETERS! ###
echo "Save Parameters in ad hoc file"
echo -e "wd=${wd}\ntag=${tag}\ninput=${input}" > "${tmp}/ENQUIRE2KG_input_parameters.txt"
echo
####

mkdir -p "$tmp"

echo

# export JOBDIR=$(pwd)/mnt-neo4j_CMITDR

neo4jlib="${tmp}/var/lib/neo4j"
neo4jlog="${tmp}/var/log/neo4j"

mkdir -p "${neo4jlib}/plugins"
mkdir -p "${neo4jlib}/data"
mkdir -p "${neo4jlib}/run"

mkdir -p "$neo4jlog"
mkdir -p "$tmp/output"

echo -e "##### GRAPH DATABASE CONSTRUCTION STARTS #####"

apptainer run \
  --mount dst=/var/lib/neo4j,src="$neo4jlib" \
  --mount dst=/var/log/neo4j,src="$neo4jlog" \
  "$sif" make_graphDB_V4.py --inputfolder "$input" --db-name neo4j --outputfolder "$tmp/output"

echo -e "##### GRAPH DATABASE CONSTRUCTION COMPLETED #####"

echo -e "##### START GRAPH DATABASE LOCATED IN ${tmp} #####"

apptainer exec \
  --mount dst=/var/lib/neo4j,src="$neo4jlib" \
  --mount dst=/var/log/neo4j,src="$neo4jlog" \
  "$sif" neo4j console