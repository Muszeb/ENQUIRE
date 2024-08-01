# ENQUIRE

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10692274.svg)](https://doi.org/10.5281/zenodo.10692274)

The accelerating growth of scientific literature overwhelms our capacity to manually distil complex phenomena like molecular networks linked to diseases. Moreover, biases in biomedical research and database annotation limit our interpretation of facts and generation of hypotheses. ENQUIRE (Expanding Networks by Querying Unexpectedly Inter-Related Entities) offers a time- and resource-efficient alternative to manual literature curation and database mining. ENQUIRE reconstructs and expands co-occurrence networks of genes and biomedical ontologies from user-selected input corpora and network-inferred PubMed queries. The integration of text mining, automatic querying, and network-based statistics mitigating literature biases makes ENQUIRE unique in its broad-scope applications. For example, ENQUIRE can generate co-occurrence gene networks that reflect high-confidence, functional networks. When tested on case studies spanning cancer, cell differentiation and immunity, ENQUIRE identified interlinked genes and enriched pathways unique to each topic, thereby preserving their underlying diversity. ENQUIRE supports biomedical researchers by easing literature annotation, boosting hypothesis formulation, and facilitating the identification of molecular targets for subsequent experimentation.

![](https://github.com/Muszeb/ENQUIRE/blob/main/ENQUIRE_graphical_abstract.png)

- If you find ENQUIRE useful to pursue your research, [please cite us](https://www.biorxiv.org/content/10.1101/2023.09.10.556351v1)

<details><summary>INSTALLATION</summary> 

ENQUIRE can currently be run on LINUX systems and LINUX virtual machines using [Apptainer/Singularity](https://apptainer.org/docs/user/latest/introduction.html). Please follow the installation steps for [Linux](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages) or [Windows/Mac](https://apptainer.org/docs/admin/main/installation.html#installation-on-windows-or-mac) install Apptainer/Singularity in order to use ENQUIRE. The file called `ENQUIRE.sif` (1.5 GB in size) is a compressed Singularity Image File (SIF) that already contains all the code, dependendencies and stable metadata needed to run ENQUIRE, so no extra installation steps are needed. We recommend adding the path to the `apptainer` executable to your `PATH` variable (e.g. by editing your `.bashrc` file). This allows to directly execute `ENQUIRE.sif` as any other executable (`./ENQUIRE.sif`).

Next, clone the repository:

```bash
git clone https://github.com/Muszeb/ENQUIRE.git
cd ENQUIRE
```

then, download the SIF image file `ENQUIRE.sif` from [FigShare](https://figshare.com/articles/software/ENQUIRE/24434845) and place it in the repository, check that the file is intact with `md5sum`, and make it executable

```
md5sum -c md5sum_ENQUIRE_sif.txt
chmod +x ENQUIRE.sif
```

You can then place the `ENQUIRE` directory or `ENQUIRE.sif` wherever you wish to, and possibly add its location to your `PATH` variable for an easier calling.

</details>

<details><summary>USAGE</summary> 

**The exemplary code snippets assume that `apptainer` location is added to your `PATH` variable, and that you're running the commands from the `ENQUIRE` main directory (do `cd /path/to/ENQUIRE` to test them).**
Here is how you call ENQUIRE scripts using `ENQUIRE.sif`:

```bash
# assuming the `apptainer` location is in your PATH variable and you did `cd ENQUIRE` or `ENQUIRE.sif` is in your working directory
Usage: ./ENQUIRE.sif <script_name> [script_argument]
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
    
- Alternatively, we also offer a Python script to extract the PubMed identifiers of all papers cited in a reading of interest (e.g. a review paper of a particular topic). From the `ENQUIRE` folder and virtual environment, type on the command line:

```bash
# assuming the `apptainer` location is in your PATH variable and you did `cd ENQUIRE` or `ENQUIRE.sif` is in your working directory
./ENQUIRE.sif efetch_references.py tag ref1 ref2 ref3 ...
```
where `tag` is the name of the plain text output file, while `ref1 ref2 ref3 ...` are the PMIDs of the papers you want to extract the references from. The output will look like the example from the previous section and is therefore ready to be used as ENQUIRE input. 
**DISCLAIMER**: if the references are not annotated into the Pubmed's API, The script will silently return no match - this may go unnoticed when fetching references from multiple articles. As a rule of thumb, look for "References" in the "page navigation" menu on the Pubmed page of the article of interest to tell the web-annotation status of an article.

</details>

<details><summary>LAUNCHING ENQUIRE</summary> 

- Before running an actual task, take a look at `ENQUIRE_methods_overview.png`: the figure briefly illustrates the main steps of the algorithm.

- In the next exemplary code snippet, we assumed you cloned this repository and `ENQUIRE` is your current working directory.

- **IMPORTANT NOTE**: it is highly recommended to get an **NCBI API_KEY** before running ENQUIRE. [Getting one is very easy](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us). You can then copy the API key and enter it as an environmental variable on the command line, like so: 
```bash
export NCBI_API_KEY=your_api_key_here
```
This will ensure your API KEY is passed as an environmental variable to all ENQUIRE runs within the same terminal session. 

- you can inspect the code Help section by running (from the `ENQUIRE` directory) `./ENQUIRE.sif ENQUIRE.sh -h`:
 
```
####################################################################################

Expanding Networks by Querying Unexpectedly Inter-Related Entities

####################################################################################

####################################################################################

Usage: ./ENQUIRE.sif ENQUIRE.sh [script_arguments]

Legend:	[-flag_short|--flag_long|config file variable, if available]:

[-p|--path|wd] = the path to the working directory (wd), where the output directory will be written in.
	It must be the ENQUIRE main folder, with ./code and ./input as subfolders.
	The default is the current working directory.

[-i|--input|to_py] = input.txt: a 'seed' input text file containing one PMID per line.
	It can be obtained from a PubMed querying specifying 'PMID' as the download format option.
	A minimun of 3 entries is required, but a list at least a few dozens articles is highly recommended.

[-t|--tag|tag] = A tag definining the task.
	It must be an alphanumeric string (underline_spaced_words are accepted).

[-j|--ncores|ncores] = The max number of CPU cores to be used.
	Default is 6.

[-c|--combine-set|comb] = how many k entities should be intersected to construct a query?
	3: loose searches, 4: moderate (default), 5: very strict queries.

[-r|--representativeness|thr] = representativeness threshold (%) for a subgraph to be included in the network expansion steps (default: 1 %).
	Example: if a subgraph contains nodes exclusively mentioned in 1 paper out of a total of 100, that subgraph has a 1% representativeness.

[-a|--attempts|A] = how many query attempts (i.e. k-sized graphlets) should be run to connect any two network communities?
	1: conservative, 2: moderate (default), 3: greedy.

[-k|--connectivity|K] = minimal community connectivity (K), which applies to any expansion-derived entities:
	each gene/MeSH term must be connected to at least K original communities to be incorporated in the expanded network - default: 2.

[-e|--entity|etype] = which entity type ('gene','MeSH') are you interested into? Omit or 'all' to textmine both entities.

[-f|--config] = if a config file is being used, specify its full path (e.g. input/textmining_config.txt).
	This option overwrites any parameter set by a different option.

[-w|--rscript|rscript] = path to the Rscript compiler. If using 'ENQUIRE.sif', it defaults to the containerized version of R.

[-d|--inputdata|sd] = path to the input data folder. If using 'ENQUIRE.sif', it defaults to the containerized input folder.
	WARNING: this option is still under development, to allow users to set different species targets
	and subsequently change the H.s. specific metadata.

[-m|--cellentitymodule|CELLTAGSBOOL] = Boolean, enable removing of character spans tagged as cell lines or types (e.g. 'CD8+ T-cell')?
	Default: False.

[-h|--help] = print this help message.

You might be seeing this Help because of an input error.

####################################################################################
``` 

Let's set up an example: we want to extract biomedical information from publications dealing with chemically-induced colitis in melanoma patients undergoing checkpoint-inhibitors therapy. Our ENQUIRE job might then look something like

```bash
# assuming the `apptainer` location is in your PATH variable and you did `cd ENQUIRE` or `ENQUIRE.sif` is in your working directory
./ENQUIRE.sif ENQUIRE.sh -t ICI_and_Colitis -i test_input/pmid-ICI_and_Colitis.txt
```

Where all the other parameters described in the `Help` message of `ENQUIRE.sh` are set to default values. The passing of the parameters could be easen by using the `ENQUIRE_config.txt` file that resides in the main `ENQUIRE` directory: the left hand side of each variable assignment must be kept unchanged, while the right hand side can be tweaked according to one's needs. Additional information on the parameters are given in `ENQUIRE_flowchart.png`. Then, the program can be launched by running:

```bash
# assuming the `apptainer` location is in your PATH variable and you did `cd ENQUIRE` or `ENQUIRE.sif` is in your working directory
./ENQUIRE.sif ENQUIRE.sh -f ENQUIRE_config.txt
```
</details>

<details><summary>EXPLANATION OF THE OUTPUT DATA STRUCTURE</summary>

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

</details>

<details><summary>EXECUTING POST-HOC ANALYSES</summary> 

#### Context-aware gene set annotation 
- Run `./ENQUIRE.sif context_aware_gene_sets.R [options]` to perform automatic annotation of gene sets, using ENQUIRE-generated, Gene/MeSH edge and node tables and Fuzzy-C-Means (FCM). See the original manuscript for further information.

```
Usage: ./ENQUIRE.sif context_aware_gene_sets.R [options]

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

	-r PARAMETER, --round=PARAMETER
		Should membership degrees be rounded to the first significant digit (helps the stability of the results)?
		default: True [T,F]

	-s PARAMETER, --setsize=PARAMETER
		minimal gene set size (default: 2)

	-v VARIANCE, --varthreshold=VARIANCE
		Dimensionality reduction based on the chosen proportion of Variance
		observed upon PCA-transforming the inverse-log-similarity between nodes (default: 0.99. range [0-1]).
		Set it to 1 to use untrasformed, scaled node similarities.

	-m MESH, --meshxgs=MESH
		How many MeSH terms which are closest to the cluster centroids should be used to describe a gene set? (default:3)

	-p PATH, --netpathdata=PATH
		Path to 'ENQUIRE-KNet_STRING_RefNet_Reactome_Paths.RData.gz' (required).
		If using the ENQUIRE.sif singularity image, the default path should point to the containerized copy of the file.

	-h, --help
		Show this help message and exit
```

- You can use the exemplary output files contained in `tmp-Ferroptosis_and_Immune_System` to test the script:
```bash
# assuming the `apptainer` location is in your PATH variable and you did `cd ENQUIRE` or `ENQUIRE.sif` is in your working directory
./ENQUIRE.sif context_aware_gene_sets.R -e tmp-Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System_Complete_edges_table_subgraph.tsv 
-n tmp-Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System_Complete_nodes_table_subgraph.tsv
```
The output will be saved in the default-tagged spreadsheet file `ENQUIRE_context_aware_gene_sets.xlsx` as well as a plot showing the reconstructed gene sets as a PNG image. Please note that the script might last quite long, due to the FCM algorithm.

#### Context-aware pathway enrichment analysis
- Run `./ENQUIRE.sif context_aware_pathway_enrichment.R [options]` to perform topology-based, pathway enrichment analysis using [SANTA](https://www.bioconductor.org/packages/devel/bioc/vignettes/SANTA/inst/doc/SANTA-vignette.html), Reactome *H. sapiens* pathways, and STRING's *H. sapiens*, physical PPI network, using ENQUIRE-generated, gene-gene edge table. See the original manuscript for further information.

```
Usage: Rscript code/context_aware_pathway_enrichment.R [options]

Options:
	-w PATH, --directory=PATH
		Working directory (default to current working directory)

	-o PATH, --outdirectory=PATH
		Output directory (default to current working directory, and must preexist)

	-n PATH, --netpathdata=PATH
		Path to 'ENQUIRE-KNet_STRING_RefNet_Reactome_Paths.RData.gz' (required).
		If using the ENQUIRE.sif singularity image, the default path should point to the containerized copy of the file.

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
		Default: holm.

	-q QSCORENET, --qscorenet=QSCORENET
		Do you want to save a copy of the STRING network in GRAPHML format with ENQUIRE-inferred QScores as node weights?
		default: False [T,F]

	-h, --help
		Show this help message and exit
```

- You can use the exemplary output files contained in `tmp-Ferroptosis_and_Immune_System` to test the script (we reduce the number of tested pathways with the `s` parameter to speed up the process):
```bash
# assuming the `apptainer` location is in your PATH variable and you did `cd ENQUIRE` or `ENQUIRE.sif` is in your working directory
./ENQUIRE.sif context_aware_pathway_enrichment.R -e tmp-Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System_Genes_edges_table_subgraph.tsv -s 30
```
The output will be saved in the default-tagged spreadsheet file `ENQUIRE_context_aware_pathway_enrichment.xlsx`, together with two PNG images showing the test statistics p-value distribution and the correlation between the Node score and degree. Please note that the script might take quite long to finish, and it benefits from a high performance computer, if available. 

</details>

<details><summary> POSSIBLE SOURCES OF ERRORS </summary>

- Test the command `which apptainer`: if `apptainer` location is not in your `PATH` variable, you need to invoke it by specifying its path, that is doing `/path/to/apptainer run /path/to/ENQUIRE.sif ...` instead of `./path/to/ENQUIRE.sif ...`

- Test the command `awk '/MemAvailable/ {print $2}' /proc/meminfo` on your command line: this is the way ENQUIRE checks the available RAM on Linux systems, in order to avoid overflows.  Make sure `awk` is installed on your system. If you witness a non-awk related issue, contact us with information on your system and possible solutions to alternatively track the available memory on your OS.

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

</details>

<details><summary> REPRODUCIBILITY </summary>

Two identical runs of ENQUIRE should produce identical co-occurrence networks and query formulations, as long as NCBI made no updates on the MeSH indexing of PubMed articles involved during the time that separates the two runs. In that case, the later run should produce queries that are supersets of the earlier one.
The exemplary output directory `tmp-Ferroptosis_and_Immune_System` was generated between 10.10.23 and 11.10.23 and has been used to generate the results illustrated in the ENQUIRE manuscript. The output was found to be reproducible on 3 different Linux Machines (2 Ubuntu and 1 ARCH-LINUX distributions).
The use of a containerized image (the SIF file) should guarantee the reproducibility irrespective of the host operating system. While several other tests on different operating systems show consistency in the network reconstruction steps, we cannot rule out the possibility that the network expansion step might diverge in some cases, irrespective of the internally coded, fixed seeds.

</details>

<details><summary>IMPORTANT INFORMATION ON PUBMED ACCESSIBILITY</summary>
	
As of 21.11.22, [important changes](https://www.nlm.nih.gov/pubs/techbull/so22/so22_updated_pubmed_e_utilities.html) have been applied to NCBI's e-utilities. In particular, it is now impossible to stream all records exceeding 10,000 PMIDs from any particular query to the PubMed database. This required to redesign the use of the e-utilities. While it's overall functionality was still preserved, we cannot guarantee the retrieval of all matching records, if the network-based queries obtained by intersecting relevant entities match more than 10,000 records (typically, this is a rare event when intersecting at least 4 distinct entities).
</details>

<details><summary>TESTED OPERATING SYSTEMS </summary>

 Below is a list of operating systems tested for installation and running of Singularity/Apptainer and ENQUIRE:

 - Linux 6.4.12-arch1-1 #1 SMP PREEMPT_DYNAMIC (x86_64 GNU/LINUX)
 - Linux 5.15.0-84-generic #93~20.04.1-Ubuntu SMP (x86_64 GNU/LINUX)
 - Virtual Machine created using Oracle Virtual Box and running Ubuntu 20 LTS 
 
</details>
