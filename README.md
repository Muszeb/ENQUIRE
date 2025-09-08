![](https://github.com/Muszeb/ENQUIRE/blob/ENQUIRE-Docker/ENQUIRE_2025_LOGO_github.png)

# ENQUIRE [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12734778.svg)](https://doi.org/10.5281/zenodo.12734778) ![GitHub License](https://img.shields.io/github/license/Muszeb/ENQUIRE)

[//]: # "<table>"
[//]: # "<tr><th> References </th><th> Distribution </th></tr>"
[//]: # "<tr><td>"

| References | Link |
| :---: | :---: |
| Methods | [<img src="https://upload.wikimedia.org/wikipedia/commons/7/77/Open_Access_logo_PLoS_transparent.svg" width="20"/>](https://doi.org/10.1371/journal.pcbi.1012745) | 
| Application | [<img src="https://upload.wikimedia.org/wikipedia/commons/7/77/Open_Access_logo_PLoS_transparent.svg" width="20"/>](https://doi.org/10.1038/s41598-025-11944-5) |  
| Source Code | [<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/d/df/Figshare_logo.svg/2560px-Figshare_logo.svg.png" width="100"/>](https://doi.org/10.1371/journal.pcbi.1012745.s001) | 
| Updates (Latest Release) | [![GitHub Release](https://img.shields.io/github/v/release/muszeb/enquire?label=%20)](https://github.com/Muszeb/ENQUIRE/releases/latest) |
| Media | [![Static Badge](https://img.shields.io/badge/https%3A%2F%2Fimg.shields.io%2Fbadge%2F-Flash_Talk-white?style=plastic&logo=youtube&logoColor=red&logoSize=auto&label=%20&labelColor=white&color=white)](https://www.youtube.com/watch?v=APwoza1JZNY) for ISMB/ECCB 2025 in Liverpool, UK |
| How to use | [Start here for Docker Implementation](#instruction-manual-docker-version) |

[//]: # "</td><td>"

| Implementation (<span id="Link">Link</span>) | Requires | Containerizes |
| :---: | :---: | :---: |
| [![Static Badge](https://img.shields.io/badge/https%3A%2F%2Fimg.shields.io%2Fbadge%2FDocker-latest-white?logo=Docker&logoColor=FFFFFF&label=Docker&labelColor=42a4f5)](https://hub.docker.com/r/muszeb/enquire) | ![Docker Image Size](https://img.shields.io/docker/image-size/muszeb/enquire) | <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/6/66/Openlogo-debianV2.svg/1200px-Openlogo-debianV2.svg.png" width="19"/> <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/4/4b/Bash_Logo_Colored.svg/512px-Bash_Logo_Colored.svg.png?20180723054350" width="21"/> <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/c3/Python-logo-notext.svg/1200px-Python-logo-notext.svg.png" width="20"/> <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/1b/R_logo.svg/724px-R_logo.svg.png?20240131042527" width="23"/> ![Static Badge](https://img.shields.io/badge/_-_-green?style=plastic&logo=Neo4j&logoColor=white&logoSize=auto) |
| [![Static Badge](https://img.shields.io/badge/https%3A%2F%2Fimg.shields.io%2Fbadge%2FApptainer-latest-white?logo=Figshare&logoColor=FFFFFF&label=Apptainer&labelColor=CC7700&color=FFFFFF)](https://doi.org/10.6084/m9.figshare.29357207.v2) | ![Static Badge](https://img.shields.io/badge/image_size-2_GiB-orange) | <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/6/66/Openlogo-debianV2.svg/1200px-Openlogo-debianV2.svg.png" width="19"/> <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/4/4b/Bash_Logo_Colored.svg/512px-Bash_Logo_Colored.svg.png?20180723054350" width="21"/> <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/c3/Python-logo-notext.svg/1200px-Python-logo-notext.svg.png" width="20"/> <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/1b/R_logo.svg/724px-R_logo.svg.png?20240131042527" width="23"/> ![Static Badge](https://img.shields.io/badge/_-_-green?style=plastic&logo=Neo4j&logoColor=white&logoSize=auto) |
| [![Static Badge](https://img.shields.io/badge/https%3A%2F%2Fimg.shields.io%2Fbadge%2FApptainer-original-white?logo=Figshare&logoColor=FFFFFF&label=Apptainer&labelColor=CC7700&color=FFFFFF)](https://doi.org/10.6084/m9.figshare.24434845.v10) | ![Static Badge](https://img.shields.io/badge/image_size-1.4_GiB-orange) | <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/6/66/Openlogo-debianV2.svg/1200px-Openlogo-debianV2.svg.png" width="19"/> <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/4/4b/Bash_Logo_Colored.svg/512px-Bash_Logo_Colored.svg.png?20180723054350" width="21"/> <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/c3/Python-logo-notext.svg/1200px-Python-logo-notext.svg.png" width="20"/> <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/1b/R_logo.svg/724px-R_logo.svg.png?20240131042527" width="23"/> |

## INSTRUCTION MANUAL (DOCKER VERSION)

<details><summary>INSTALLATION</summary> 

ENQUIRE can currently be run on Linux, Windows and macOS systems using [Docker](https://docs.docker.com/guides/docker-overview/). Please [install Docker](https://docs.docker.com/engine/install/) following the specific steps for your OS in order to use ENQUIRE. The Docker image necessary to run ENQUIRE is hosted on [Docker Hub](https://hub.docker.com/r/muszeb/enquire). The file contains all the code, dependendencies and stable metadata needed to run ENQUIRE, so no extra installation steps are needed.

In order to download and start using ENQUIRE, open your command line and run:
```
docker pull muszeb/enquire
```

This will pull ENQUIRE from Docker Hub. Wait until the process is finished, and check the existence of the Docker image with:
```bash
# assuming you have downloaded the ENQUIRE image like above:
docker images
```

You should see ENQUIRE listed among the images. The `REPOSITORY` and `TAG` columns will inform you about the image name and the current version of ENQUIRE.

Next, choose a directory and clone the ENQUIRE repository with Git or download the data as a zip file:

```bash
git clone https://github.com/Muszeb/ENQUIRE.git
```

Now, navigate to the `ENQUIRE` folder you just downloaded and create a Docker container of ENQUIRE giving it a name of your choice with:
```bash
# Linux and Mac; assuming you are working from the ENQUIRE folder
docker run  -v $(pwd):/mnt/ --name your_container_name -it muszeb/enquire:latest
```

```bash
# Windows; assuming you are working from the ENQUIRE folder
docker run -v %cd%:/mnt/ --name your_container_name -it muszeb/enquire:latest
```

Where `-v %cd%:/mnt/` mounts the current directory (the contents of the `ENQUIRE` folder) to the Docker container just created, such that the program can access, edit and write files in your local host.

**_NOTE:_**  Take into account that whenever we want to indicate the use of a local directory or file from inside ENQUIRE's Docker container, we must use the linked/mounted directory `/mnt`.

**_NOTE 2:_**  The Docker daemon runs as the root user by default, meaning that files ENQUIRE produces in the mounted directory `/mnt` will be written as the root user. This means that permission conflicts might likely arise whenever ENQUIRE is run in a cluster with mutiple users. In order for a user to use Docker, they will necessarily receive root-level privileges, which will impact the security of your multi-user system. The proper handling of these limitations on Docker's side is beyond the scope of this README. For more information, please refer to the `Multi-user systems` section of this README.

You should see a terminal open up with `/mnt/` as a working directory. You may check the mounting worked well with:
```bash
# check files were correctly linked inside the Docker container
ls
```

You should see the files inside of the `ENQUIRE` folder.

You can exit ENQUIRE's Docker container any time with `Control + D` or by typing `exit` on the command line. To start ENQUIRE again, you can start the already existing Docker container with:
```bash
# assuming you created a container with the name your_container_name
docker start -i your_container_name
```

If you do not remember the name given to the container, you can always check with:
```bash
# assuming you created a container with the name your_container_name
docker ps -all
```

You are then ready to use ENQUIRE.

[Back to the beginning of the instruction manual](#instruction-manual-docker-version)

</details>

<details><summary>USAGE</summary> 

**The exemplary code snippets assume that you're that you have already loaded the ENQUIRE Docker image and that you are inside the bash terminal inside of the Docker container.**

Here is how you call ENQUIRE scripts using Docker:

```bash
# assuming you have already loaded the ENQUIRE Docker image, and that you are running the commands from the ENQUIRE main directory
Usage: run.sh <script_name> [script_argument]
```

Where `<script_name>` is one of:

- `efetch_references.py`
- `ENQUIRE.sh`
- `context_aware_gene_sets.R`
- `context_aware_pathway_enrichment.R`

</details>

<details><summary>INPUT FILE</summary>

A valid input file should consist of a list of PubMed Identifiers (PMIDs) stored in plain text files, one PMID per lines, such as:

 
    26250731
    22835603   
    31254155
    32658557
    30729513
    31620854
    30338457
    33711241
    28640701
    24725689

- The easiest way to generate a valid ENQUIRE input file is to generate a [PubMed query on the NCBI's website](https://pubmed.ncbi.nlm.nih.gov/). Use of MeSH terms and exclusion of review articles is recommended but not mandatory. Then, click on **Save**, choose **Selection: All results** and **Format: PMID**, and **Create file**: 
![Exemplary PubMed Query with ENQUIRE-compliant Save options](https://github.com/Muszeb/ENQUIRE/blob/main/Example_Input_PubMed_Query.png)
    
- Alternatively, we also offer a Python script to extract the PubMed identifiers of all papers cited in a reading of interest (e.g. a review paper of a particular topic). From the `ENQUIRE` Docker container, type on the command line:

```bash
run.sh efetch_references.py tag ref1 ref2 ref3 ...
```
where `tag` is the name of the plain text output file, while `ref1 ref2 ref3 ...` are the PMIDs of the papers you want to extract the references from. The output will look like the example from the previous section and is therefore ready to be used as ENQUIRE input. 
**DISCLAIMER**: if the references are not annotated into the Pubmed's API, The script will silently return no match - this may go unnoticed when fetching references from multiple articles. As a rule of thumb, look for "References" in the "page navigation" menu on the Pubmed page of the article of interest to tell the web-annotation status of an article.

[Back to the beginning of the instruction manual](#instruction-manual-docker-version)

</details>

<details><summary>LAUNCHING ENQUIRE</summary> 

- Before running an actual task, take a look at `ENQUIRE_methods_overview.png`: the figure briefly illustrates the main steps of the algorithm.

- In the next exemplary code snippet, we assume you have followed the installation steps and are working from inside of the bash terminal in ENQUIRE's Docker container.

- **IMPORTANT NOTE**: it is highly recommended to get an **NCBI API_KEY** before running ENQUIRE. [Getting one is very easy](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us). You can then copy the API key and enter it as an environmental variable on the command line, like so: 
```bash
export NCBI_API_KEY=your_api_key_here
```
This will ensure your API KEY is passed as an environmental variable to all ENQUIRE runs within the same terminal session. 

- you can inspect the code Help section by running (from the `ENQUIRE` directory) `run.sh ENQUIRE.sh -h`:
 
```
####################################################################################

Expanding Networks by Querying Unexpectedly Inter-Related Entities

####################################################################################

####################################################################################

Usage: run.sh ENQUIRE.sh [script_arguments]

Legend:	[-flag_short|--flag_long|config file variable, if available]:

[-p|--path|wd] = the path to the working directory (wd), where the output directory will be written in.
	It must be the ENQUIRE main folder, with ./code and ./input as subfolders.
	The default is the current working directory.

[-i|--input|to_py] = input.txt: a 'seed' input text file containing one PMID per line.
	It can be obtained from a PubMed querying specifying 'PMID' as the download format option.
	A minimun of 3 entries is required, but a list at least a few dozens articles is highly recommended.

[-t|--tag|tag] = A characteristic tag definining the task.
	It must be an alphanumeric string.

[-j|--ncores|ncores] = The max number of CPU cores to be used.
	Default is 6.

[-c|--combine-set|comb] = how many N entities to intersect to construct a query?
	3: loose searches, 4: moderate (default), 5: very strict queries.

[-r|--representativeness|thr] = representativeness threshold (%) for a subgraph to be included in the network expansion steps? (default: 0 %).
	Example: if a subgraph contains nodes exclusively mentioned in 10 papers out of a total of 100, that subgraph has a 10% representativeness.

[-a|--attempts|A] = how many query attempts (i.e. pairs of motifs or genes) should be run in order to connect any two subgraphs?
	1: conservative, 2: moderate (default), 3: greedy.

[-k|--connectivity|K] = minimal community connectivity (K), which applies to any expansion-derived entities:
	each gene/MeSH term must be connected to at least K original communities to be incorporated in the expanded network - default: 2.

[-e|--entity|etype] = which entity type (gene/MeSH) are you interested into? Omit or 'all' to textmine both entities.

[-f|--config] = if a config file is being used, specify its full path (e.g. input/textmining_config.txt).
	This option overwrites any parameter set by a different option.

[-w|--rscript|rscript] = path to the Rscript compiler. It defaults to the containerized version of R.

[-d|--inputdata|sd] = path to the input data folder compiler. It defaults to the containerized input folder.
	WARNING: this option is still under development, to allow users to set different species targets
	and subsequently change the H.s. specific metadata.

[-m|--cellentitymodule|CELLTAGSBOOL] = Boolean, enable removing of character spans tagged as cell lines or types (e.g. 'CD8+ T-cell')?
	Default: False.

[-h|--help] = print this help message.

You might be seeing this Help because of an input error.

####################################################################################
``` 

Let's set up an example: we want to know the current state-of-the-art regarding chemically-induced colitis in melanoma patients undergoing checkpoint-inhibitors therapy. Our ENQUIRE job might then look something like

```bash
run.sh ENQUIRE.sh -t ICI_and_Colitis -i /mnt/test_input/pmid-ICI_and_Colitis.txt
```

Where all the other parameters described in the `Help` message of `ENQUIRE.sh` are set to default values. The passing of the parameters could be easen by using the `ENQUIRE_config.txt` file that resides in the main `ENQUIRE` directory: the left hand side of each variable assignment must be kept unchanged, while the right hand side can be tweaked according to one's needs. Additional information on the parameters are given in `ENQUIRE_flowchart.png`. Then, the program can be launched by running:

```bash
# assuming you are working from within ENQUIRE's Docker container.
run.sh ENQUIRE.sh -f ENQUIRE_config.txt
```
[Back to the beginning of the instruction manual](#instruction-manual-docker-version)

</details>

<details><summary>EXPLANATION OF THE OUTPUT DATA STRUCTURE</summary>

**_NOTE:_**  The Docker daemon runs as the root user by default, meaning that files ENQUIRE produces in the mounted directory `/mnt` will be written as the root user. If you are having trouble opening the files due to permission issues, please refer to the `Multi-user systems` section of this README.

- Provided a recognisable `tag` has been passed to textmining algorithm, a typical output would produce a folder named `tmp-tag`, which in turn contains as many subdirectories as the number of steps/iterations performed. For example, if the algorithm performed 
    
    1. Reconstruction of a Gene/Mesh network from the original set of papers;
    2. One query expansion and network reconstruction as the Gene/Mesh network was not fully connected yet;
    3. One query expansion and network reconstruction as the gene-gene network was not fully connected yet, then stopped;

    Then there will be three subfolders, namely `tag`, `tag_subgraph_expansion1`, `tag_subgraph_expansion2`. The counter attached to folders and file names records the subsequent attempts to the expansion and reconstruction of co-occurence networks.

	![](https://github.com/Muszeb/ENQUIRE/blob/main/output_overview/main_structure.png)

   Typically, within each of these sub-folders/iterations, three pairs of edge and node tables can be found, respectively corresponding to "Complete" (Gene/Mesh), "Gene"- and "Mesh"-only networks (TSV files). These files can be easily imported in Cytoscape or similar graph visualization tools.

	![](https://github.com/Muszeb/ENQUIRE/blob/main/output_overview/node_edgetabs.png)

	Whenever it wasn't possible to obtain one or more of the aforementioned networks, the pipeline should print a message with information on the most meaningful files to look at. It is worth mentioning that the file `tag...Complete_literature_links.tsv` within each subfolder allows fast retrieval of specific edge-associated papers by means of encoded hyperlinks.

 	![](https://github.com/Muszeb/ENQUIRE/blob/main/output_overview/litlinks.png)

	The batch of queries that were tested in each iteration is stored in `tag...ordered_queries.tsv` within each respective subfolder. Additional meta-data can be explored under the `data/` subfolder. Besides node and edge tables for individual subgraphs (i.e. gene/MeSH of gene-only connected components), here you could also explore how the original co-occurrence multigraph looked like, before the network-based test statistics (`tag...edge_list_allxall.tsv`). 
	
	![](https://github.com/Muszeb/ENQUIRE/blob/main/output_overview/data.png)

    Furthemore, under `tmp-tag`, the file `source_pmids.txt` contains all the inspected articles for the given ENQUIRE job. These can also be consulted specifically for each iteration under `tmp-tag/efetch_inputs`.
   
    Please don't hesitate to contact us for any clarification on the purposes of any file.
  
- Interactive .html networks

	It is also possible to visually inspect Gene-MeSH networks and the reduced networks containing only cliques in two .html files, respectively stored within each iteration's subfolder as `tag...interactive_Gene-MeSH_Network.html` and `tag...interactive_Cliques_Network.html`.

	![](https://github.com/Muszeb/ENQUIRE/blob/main/output_overview/html.png)

[Back to the beginning of the instruction manual](#instruction-manual-docker-version)

</details>

<details><summary>EXECUTING POST-HOC ANALYSES</summary> 

#### Context-aware gene set annotation 
- Run `run.sh context_aware_gene_sets.R [options]` to perform automatic annotation of gene sets, using ENQUIRE-generated, Gene/MeSH edge and node tables and Fuzzy-C-Means (FCM). See the original manuscript for further information.

```
Usage: run.sh context_aware_gene_sets.R [options]

Options:
	-w PATH, --directory=PATH
		Output directory [default to current working directory]

	-e PATH, --edgetable=PATH
		Path to an ENQUIRE-generated, Gene/MeSH edge table file (required)

	-n PATH, --nodetable=PATH
		Path to an ENQUIRE-generated, Gene/MeSH node table file (required)

	-t TAG, --tag=TAG
		tag prefix (default to 'ENQUIRE')

	-d PARAMETER, --membdeg=PARAMETER
		minimal membership degree for gene-to-cluster association (default: 0.05), range [0-1]

	-s PARAMETER, --setsize=PARAMETER
		minimal gene set size (default: 2)

	-h, --help
		Show this help message and exit
```

- You can use the exemplary output files contained in `tmp-Ferroptosis_and_Immune_System` to test the script:
```bash
run.sh context_aware_gene_sets.R -e tmp-Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System_Complete_edges_table_subgraph.tsv -n tmp-Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System_Complete_nodes_table_subgraph.tsv
```
Please note that the script might last quite long, due to the FCM algorithm.

#### Context-aware pathway enrichment analysis
- Run `eun.sh context_aware_pathway_enrichment.R [options]` to perform topology-based, pathway enrichment analysis using [SANTA](https://www.bioconductor.org/packages/devel/bioc/vignettes/SANTA/inst/doc/SANTA-vignette.html), Reactome *H. sapiens* pathways, and STRING's *H. sapiens*, physical PPI network, using ENQUIRE-generated, gene-gene edge table. See the original manuscript for further information.

```
Usage: Rscript code/context_aware_pathway_enrichment.R [options]

Options:
	-w PATH, --directory=PATH
		Working directory (default to current working directory)

	-o PATH, --outdirectory=PATH
		Output directory (default to current working directory, and must preexist)

	-n PATH, --netpathdata=PATH
		Path to 'ENQUIRE-KNet_STRING_RefNet_Reactome_Paths.RData.gz' (required).
		If the current working directory is not the 'ENQUIRE' folder, the default path ('input/...') will throw an error.

	-e PATH, --edgetable=PATH
		Path to an ENQUIRE-generated, gene-gene edge table file (required).

	-c PARAMETER, --cores=PARAMETER
		max number of cores used (PSOCK parallelization) (default: 4), >1 recommended.

	-t TAG, --tag=TAG
		tag prefix (default to 'ENQUIRE').

	-s PARAMETER, --setsize=PARAMETER
		maximum Reactome pathway size (default: 100, minimum 3).

	-p PARAMETER, --permutations=PARAMETER
		number of permutations to infer KNet null distribution
		(default: 100, the higher the more accurate the test statistics).

	-f PARAMETER, --padjust=PARAMETER
		P-value adjustment method, must be one of [holm, hochberg, hommel, bonferroni, BH, BY, fdr, none].
		Default and recommended: holm, as the p-value null distribution is not guaranteed to be uniform.

	-h, --help
		Show this help message and exit
```

- You can use the exemplary output files contained in `tmp-Ferroptosis_and_Immune_System` to test the script:
```bash
run.sh context_aware_pathway_enrichment.R -e tmp-Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System_Genes_edges_table_subgraph.tsv -s 30
```
Please note that the script might last quite long, and it benefits from a high performance computer, if available. 

[Back to the beginning of the instruction manual](#instruction-manual-docker-version)

</details>

<details><summary> POSSIBLE SOURCES OF ERRORS </summary>

- Test the command `mawk '/MemAvailable/ {print $2}' /proc/meminfo` on the command line of the Docker container: this is the way ENQUIRE checks the available RAM on Linux systems, in order to avoid overflows. If you witness a non-awk related issue, contact us with information on your system and possible solutions to alternatively track the available memory on your OS.

- When computing large networks, an error related to the default `Stack Size` can potentially appear, especially when running R scripts, such as `Error: C stack usage is too close to the limit`. In this case, one shall set a higher stacksize to allow the script to complete, via 

    ```
    ulimit -s N 
    ```    
    Where `N` shall be a size expressed in Kb to set as the maximum stack size. You could first check the number returned by `Cstack_info()` in an active R shell. You can read more about the issue [here](https://stackoverflow.com/questions/14719349/error-c-stack-usage-is-too-close-to-the-limit) and [here](https://rdrr.io/r/base/Cstack_info.html).

- If you get a `curl`-related error of the form 
```bash
HTTP/1.1 400 Bad Request
 WARNING:  FAILURE ( Thu Feb 15 10:24:24 AM CET 2024 )
```

It means that NCBI is not willing to process your request. Sometimes, this can be due to a server hiccup, but most times using an API KEY fixes the issue. [Getting one is very easy](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us). You can then copy the API key and enter it as an environmental variable on the command line, like so: 
```bash
export NCBI_API_KEY=your_api_key_here
```
This will ensure your API KEY is passed as an environmental variable to all ENQUIRE runs within the same terminal session. 

[Back to the beginning of the instruction manual](#instruction-manual-docker-version)

</details>

<details><summary> REPRODUCIBILITY </summary>

Two identical runs of ENQUIRE should produce identical co-occurrence networks and query formulations, as long as NCBI made no updates on the MeSH indexing of PubMed articles involved during the time that separates the two runs. In that case, the later run should produce queries that are supersets of the earlier one.
The exemplary output directory `tmp-Ferroptosis_and_Immune_System` was generated between 10.10.23 and 11.10.23 and has been used to generate the results illustrated in the ENQUIRE manuscript. The output was found to be reproducible on 3 different Linux Machines (2 Ubuntu and 1 ARCH-Linux distributions).
The use of a containerized image should guarantee the reproducibility irrespective of the host operating system. While several other tests on different operating systems show consistency in the network reconstruction steps, we cannot rule out the possibility that the network expansion step might diverge in some cases, irrespective of the internally coded, fixed seeds.

[Back to the beginning of the instruction manual](#instruction-manual-docker-version)

</details>

<details><summary>MULTI-USER SYSTEMS</summary> 
As mentioned in the notes in the `installation` section and in the `explanation of the output` section, the Docker daemon runs as the root user by default, meaning that any user using ENQUIRE in a multi-user system will inevitably have access to root/admin level privileges. This impacts the security of your multi-user system. For more information on this, you can read:

[Manage Docker as a non-root user](https://docs.docker.com/engine/install/Linux-postinstall/#manage-docker-as-a-non-root-user).

To actually allow the Docker Daemon to run rootless to address these security concerns, please refer to:
[Run the Docker daemon as a non-root user](https://docs.docker.com/engine/security/rootless/).

As an alternative to Docker, you may also run images with [Podman](https://podman.io/), which also allows the setup of multiple users:
[Basic Setup and Use of Podman in a Rootless environment](https://github.com/containers/podman/blob/main/docs/tutorials/rootless_tutorial.md).

If you decide to run ENQUIRE with a root-based Docker daemon, the files outputted by ENQUIRE will be under root/admin permissions. This might prevent you from accessing files you have produced when you exit ENQUIRE's terminal. In order to address this issue, we offer a quick fix. Find your username and group name identifier by running:
```bash
# For username, groupname in Linux and MacOS.
id -u # username
id -g # groupname
```

```bash
# For Windows.
whoami # username
whoami /groups # groups the user is part of
```

Once you have retrieved your username and groupname identifiers, you may grant yourself permission

```bash
# To be executed inside ENQUIRE's terminal (which is run as root)
chown -R username:groupname /mnt
```
[Back to the beginning of the instruction manual](#instruction-manual-docker-version)

</details>

<details><summary>IMPORTANT INFORMATION ON PUBMED ACCESSIBILITY</summary>
	
As of 21.11.22, [important changes](https://www.nlm.nih.gov/pubs/techbull/so22/so22_updated_pubmed_e_utilities.html) have been applied to NCBI's e-utilities. In particular, it is now impossible to stream all records exceeding 10,000 PMIDs from any particular query to the PubMed database. This required to redesign the use of the e-utilities. While it's overall functionality was still preserved, we cannot guarantee the retrieval of all matching records, if the network-based queries obtained by intersecting relevant entities match more than 10,000 records (typically, this is a rare event when intersecting at least 4 distinct entities).

[Back to the beginning of the instruction manual](#instruction-manual-docker-version)

</details>

<details><summary>TESTED OPERATING SYSTEMS </summary>

 Below is a list of operating systems tested for installation and running of the dockerized version of ENQUIRE:

 - Linux 6.4.12-arch1-1 #1 SMP PREEMPT_DYNAMIC (x86_64 GNU/Linux)
 - Linux 5.15.0-84-generic #93~20.04.1-Ubuntu SMP (x86_64 GNU/Linux)
 - Virtual Machine created using Oracle Virtual Box and running Ubuntu 20 LTS
 - MacOS Catalina 15.7 (Docker implementation, mid-2012 MacBook Pro)
 - Windows 10 (Docker implementation) 

[Back to the beginning of the instruction manual](#instruction-manual-docker-version)
 
</details>
