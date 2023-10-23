# ENQUIRE

The accelerating growth in scientific literature is overwhelming our capacity to manually distil complex phenomena like molecular networks linked to diseases. Moreover, confirmation biases in search engines and databases influence the interpretation of facts and the generation of hypotheses. ENQUIRE (Expanding Networks by Querying Unexpectedly Inter-Related Entities) offers an alternative to manual literature curation and database mining to study complex biomedical phenomena. ENQUIRE generates a co-occurrence network of genes and biomedical ontologies (MeSH) using a corpus of publications as input. The algorithm iteratively reconstructs and expands the network by generating PubMed queries from significant interrelations, until it obtains an interconnected network of genes and MeSH relevant to the input corpus. We systematically evaluated ENQUIRE’s versatility and found that it generates co-occurrence networks similar to those based on co-expression data and manually annotated databases with a fraction of the time and human resources. Using case studies spanning cancer, cell differentiation and immunity, ENQUIRE proved to identify interlinked genes and enriched pathways unique to each topic, thereby preserving their underlying diversity. ENQUIRE supports biomedical researchers by easing literature annotation, boosting hypothesis formulation and facilitating the identification of molecular targets for subsequent experimentation.

![](https://github.com/Muszeb/ENQUIRE/blob/main/ENQUIRE_graphical_abstract.png)

- If you find ENQUIRE useful to pursue your research, [please cite us](https://www.biorxiv.org/content/10.1101/2023.09.10.556351v1)

<details><summary>INSTALLATION</summary> 

ENQUIRE can currently be run on LINUX systems and LINUX virtual machines using [Apptainer/Singularity](https://apptainer.org/docs/user/latest/introduction.html). Please follow the [installation steps](https://github.com/apptainer/apptainer/blob/main/INSTALL.md) specific to your setup and install Apptainer/Singularity in order to use ENQUIRE. The file called `ENQUIRE.sif` (1.3 GB in size) is a compressed Singularity Image File (SIF) that already contains all the code, dependendencies and stable metadata needed to run ENQUIRE, so no extra installation steps are needed.

Next, clone the repository, download the SIF image file `ENQUIRE.sif` from [our lab's website](http://sysbiomed-erlangen.weebly.com/) within the repository main folder, and make it executable:

```bash
git clone https://github.com/Muszeb/ENQUIRE.git
cd ENQUIRE 
curl https://jveralab.net/ENQUIRE.sif -o ENQUIRE.sif
chmod +x ENQUIRE.sif
```

You can then place the `ENQUIRE` directory of `ENQUIRE.sif` wherever you wish to, and possibly add its location to your `PATH` variable for an easier calling.

</details>

<details><summary>USAGE</summary> 

**The examplary code snippets assume that you're running the commands from the `ENQUIRE` main directory (do `cd /path/to/ENQUIRE` to test them).**
Here is how you call ENQUIRE scripts using `ENQUIRE.sif`:

```bash
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
./ENQUIRE.sif efetch_references.py tag ref1 ref2 ref3 ...
```
where `tag` is the name of the plain text output file, while `ref1 ref2 ref3 ...` are the PMIDs of the papers you want to extract the references from. The output will look like the example from the previous section and is therefore ready to be used as ENQUIRE input. 
DISCLAIMER: if the references are not annotated into the Pubmed's API, an error such as 

```python
File "code/efetch_references.py", line 28, in <module>
    refs+=refparse(p)
  File "code/efetch_references.py", line 20, in refparse
    refs=dpath.get(data,"**/Link") # list of {Id:value} dicts
  File "/home/musellla/miniconda3/envs/wokenv/lib/python3.8/site-packages/dpath/util.py", line 178, in get
    raise KeyError(glob)
KeyError: '**/Link'`
```
might occur. As a rule of thumb, look for "MeSH terms" in the "page navigation" menu on the Pubmed page of the article of interest.

</details>

<details><summary>LAUNCHING ENQUIRE</summary> 

- Before running an actual task, take a look at `ENQUIRE_methods_overview.png`: the figure briefly illustrates the main steps of the algorithm.

- After the download, you should see a folder called `ENQUIRE`: this is the main directory from which the program is supposed to be run.
   
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

	[-w|--rscript|rscript] = path to the Rscript compiler. If using 'ENQUIRE.sif', it defaults to the containerized version of R.

	[-d|--inputdata|sd] = path to the input data folder compiler. If using 'ENQUIRE.sif', it defaults to the containerized input folder.
		WARNING: this option is still under development, to allow users to set different species targets
		and subsequently change the H.s. specific metadata.

	[-h|--help] = print this help message.

	You might be seeing this Help because of an input error.

	####################################################################################
    ``` 

    Let's set up an example: we want to know the current state-of-the-art regarding chemically-induced colitis in melanoma patients undergoing checkpoint-inhibitors therapy. Our ENQUIRE job might then look something like

    ```bash
		# assuming you did `cd ENQUIRE` and `ENQUIRE.sif` resides there
    	./ENQUIRE.sif ENQUIRE.sh -t ICI_and_Colitis -i test_input/pmid-ICI_and_Colitis.txt
    ```

    Where all the other parameters described in the `Help` message of `ENQUIRE.sh` are set to default values. The passing of the parameters could be easen by using the `ENQUIRE_config.txt` file that resides in the main `ENQUIRE` directory: the left hand side of each variable assignment must be kept unchanged, while the right hand side can be tweaked according to one's needs. Additional information on the parameters are given in `ENQUIRE_flowchart.png`. Then, the program can be launched by running:

    ```bash
		# assuming you did `cd ENQUIRE` and `ENQUIRE.sif` resides there
    	./ENQUIRE.sif ENQUIRE.sh -f ENQUIRE_config.txt
    ```
</details>

<details><summary>EXPLANATION OF THE OUTPUT DATA STRUCTURE</summary>

- Provided a recognisable `tag` has been passed to textmining algorithm, a typical output would produce a folder named `tmp-tag`, which in turn contains as many subdirectories as the number of steps/iterations performed. For example, if the algorithm performed 
    
    1. Reconstruction of a Gene/Mesh network from the original set of papers;
    2. One query expansion and network reconstruction as the Gene/Mesh network was not fully connected yet;
    3. One query expansion and network reconstruction as the gene-gene network was not fully connected yet, then stopped;

    Then there will be three subfolders, namely `tag`, `tag_subgraph_expansion1`, `tag_subgraph_expansion2`. The counter attached to folders and file names records the subsequent attempts to the expansion and reconstruction of co-occurence networks.

	![](https://github.com/Muszeb/ENQUIRE/blob/ENQUIRE-MACOS/output_overview/main_structure.png)

   Typically, within each of these sub-folders/iterations, three pairs of edge and node tables can be found, respectively corresponding to "Complete" (Gene/Mesh), "Gene"- and "Mesh"-only networks (TSV files). These files can be easily imported in Cytoscape or similar graph visualization tools.

	![](https://github.com/Muszeb/ENQUIRE/blob/ENQUIRE-MACOS/output_overview/node_edgetabs.png)

	Whenever it wasn't possible to obtain one or more of the aforementioned networks, the pipeline should print a message with information on the most meaningful files to look at. It is worth mentioning that the file `tag...Complete_literature_links.tsv` within each subfolder allows fast retrieval of specific edge-associated papers by means of encoded hyperlinks.

 	![](https://github.com/Muszeb/ENQUIRE/blob/ENQUIRE-MACOS/output_overview/litlinks.png)

	The batch of queries that were tested in each iteration is stored in `tag...ordered_queries.tsv` within each respective subfolder. Additional meta-data can be explored under the `data/` subfolder. Besides node and edge tables for individual subgraphs (i.e. gene/MeSH of gene-only connected components), here you could also explore how the original co-occurrence multigraph looked like, before the network-based test statistics (`tag...edge_list_allxall.tsv`). 
	
	![](https://github.com/Muszeb/ENQUIRE/blob/ENQUIRE-MACOS/output_overview/data.png)

    Furthemore, under `tmp-tag`, the file `source_pmids.txt` contains all the inspected articles for the given ENQUIRE job. These can also be consulted specifically for each iteration under `tmp-tag/efetch_inputs`.
   
    Please don't hesitate to contact us for any clarification on the purposes of any file.
  
- Interactive .html networks

	It is also possible to visually inspect Gene-MeSH networks and the reduced networks containing only cliques in two .html files, respectively stored within each iteration's subfolder as `tag...interactive_Gene-MeSH_Network.html` and `tag...interactive_Cliques_Network.html`.

	![](https://github.com/Muszeb/ENQUIRE/blob/ENQUIRE-MACOS/output_overview/html.png)

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

	-s PARAMETER, --setsize=PARAMETER
		minimal gene set size (default: 2)

	-h, --help
		Show this help message and exit
```

- You can use the exemplary output files contained in `tmp-Ferroptosis_and_Immune_System` to test the script:
```bash
# assuming you did `cd ENQUIRE` and `ENQUIRE.sif` resides there
./ENQUIRE.sif context_aware_gene_sets.R -e tmp-Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System_Complete_edges_table_subgraph.tsv 
-n tmp-Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System_Complete_nodes_table_subgraph.tsv
```
Please note that the script might last quite long, due to the FCM algorithm.

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
# assuming you did `cd ENQUIRE` and `ENQUIRE.sif` resides there
./ENQUIRE.sif context_aware_pathway_enrichment.R -e tmp-Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System/Ferroptosis_and_Immune_System_Genes_edges_table_subgraph.tsv -s 30
```
Please note that the script might last quite long, and it benefits from a high performance computer, if available. 

</details>

<details><summary> POSSIBLE SOURCES OF ERRORS </summary>

- Test the command `awk '/MemAvailable/ {print $2}' /proc/meminfo` on your command line: this is the way ENQUIRE checks the available RAM on Linux systems, in order to avoid overflows.  Make sure `awk` is installed on your system. If you witness a non-awk related issue, contact us with information on your system and possible solutions to alternatively track the available memory on your OS.

- When computing large networks, an error related to the default `Stack Size` can potentially appear, especially when running R scripts, such as `Error: C stack usage is too close to the limit`. In this case, one shall set a higher stacksize to allow the script to complete, via 

    ```
    ulimit -s N 
    ```    
    Where `N` shall be a size expressed in Kb to set as the maximum stack size. You could first check the number returned by `Cstack_info()` in an active R shell. You can read more about the issue [here](https://stackoverflow.com/questions/14719349/error-c-stack-usage-is-too-close-to-the-limit) and [here](https://rdrr.io/r/base/Cstack_info.html). 

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
