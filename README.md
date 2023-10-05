# ENQUIRE

<details><summary>INSTALLATION</summary> 

ENQUIRE can currently be installed and run on UNIX systems such as MacOS, Linux and UNIX emulators/virtual machines for Windows such as WSL, CygWin or MSYS2 (installation of UNIX emulators is not included in this guide).

Start by cloning and accessing the repository. If using MacOS, clone the `ENQUIRE-MACOS` branch (this branch assumes you have MacOS-specific `coreutils` and `findutils` installed). 
```bash
# Linux/WSL/CygWing/MSYS2 
git clone https://github.com/Muszeb/ENQUIRE.git
cd ENQUIRE
```
```bash
# MacOS
git clone https://github.com/Muszeb/ENQUIRE.git
git checkout ENQUIRE-MACOS
cd ENQUIRE
```

ENQUIRE installation consists of 5 steps, although your system may already satisfy some of the requirements. 

1) Installing `curl` to be able to download and install a virtual environment manager and `EDirect`.
2) Downloading a virtual environment manager.
3) Creating a virtual environment and install Python and R package dependencies.
4) Installing `EDirect`.
5) Installing `Pandoc`.

Below is an exemplary recipe for the whole installation process. If you are using a Debian GNU/Linux system with `sudo` rights, you may install everything by doing `bash code/install_everything.sh`. Please check the installation specifications and change them accordingly to fit to your OS and working environment. Windows user may want to first [locate their local disk from a UNIX virtual machine](https://askubuntu.com/questions/943006/how-to-navigate-to-c-drive-in-bash-on-wsl-ubuntu).  

```bash
: '
1) Install `curl` to be able to download and install a virtual environment manager and EDirect.
`apt-get` is a Debian-specific package manager, you may want to check which package manager works best on your OS:
https://everything.curl.dev/get.
'
sudo apt-get update
sudo apt-get install -y curl
```
```bash
: '
2) Download a virtual environment manager. Here, we show how to download and install `micromamba`.
Suggested alternatives are `mamba` and `miniconda` - see https://mamba.readthedocs.io/en/latest/installation.html.
'
#Change the following variable to customize the installation path:
mambapath=$HOME/bin

# Download micromamba and export path to micromamba (possibly update .bashrc)
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj "${mambapath}/micromamba"
export PATH="${mambapath}:$PATH"
echo "export PATH=${mambapath}:$PATH" >> ${HOME}/.bashrc

# This command allows micromamba to execute shell commands
eval "$(micromamba shell hook --shell bash)"
```
```bash
: '
3) Create a virtual environment and install Python and R package dependencies.
Do `cd ENQUIRE` to run the following commands from within ENQUIRE main folder).
If using a different environment manager like mamba,
do `mamba env create -n ENQUIRE -f input/ENQUIRE.yml` followed by `mamba activate ENQUIRE`
and skip to the `pip install` command.
'
# Create an ad hoc environment.
# If you have another virtual environment manager installed (e.g miniconda),
# the environment might be installed under '$HOME/othermanager/envs'
micromamba create -n ENQUIRE 
micromamba activate ENQUIRE

# Install Python libraries `ENQUIRE` environment, from the `input/ENQUIRE.yml` file.
micromamba install -y -q -f input/ENQUIRE.yml
pip install -r input/PackagesNotFound.txt # these packages can't be installed by environment managers

# Install R libraries under the ENQUIRE environment (remember to "activate ENQUIRE"!)
$(which Rscript) code/install_R_libraries.R 

# Clean environment
micromamba clean --all --yes
```
```bash
: '
4) Install EDirect.
See here for the latest  and OS-specific installation command:
https://www.ncbi.nlm.nih.gov/books/NBK179288/ 
'
# Install EDirect under your HOME directory - manually adding the edirect path to .bash_profile keeps the latter cleaner
yes n | sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
echo "export PATH=$PATH:$HOME/edirect" >> $HOME/.bash_profile # necessary 
```
```bash
: '
5) Install Pandoc.
See here for the latest  and OS-specific installation command: https://pandoc.org/installing.html
'
# Install Pandoc 
sudo apt-get install -y pandoc
```

Alternatively (but not recommended), you can install all the required R and Python libraries independent of a virtual environment manager (steps 2-3). First, install `python>=3.8`, `R>=4.2`, and `pip >= 20.2`. Then, from the `ENQUIRE` directory run

```bash
pip install -r input/ENQUIRE_pip_requirements.txt
Rscript code/install_R_libraries.R
```
Note that this does not replace steps 1, 4 and 5 described above.
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
python code/efetch_references.py tag ref1 ref2 ref3 ...
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

- After the download, you should see a folder called `ENQUIRE`: this is the main directory from which the program is supposed to be run.

- Before running an actual task, take a look at `ENQUIRE_methods_overview.png`: the figure briefly illustrates the main steps of the algorithm.
   
- Now activate the `ENQUIRE` environment (unless you followed the pip-only installation) and inspect the code Help section by running (from the `ENQUIRE` directory):
 
    ```bash
    conda activate ENQUIRE
    cd ENQUIRE
    ./code/ENQUIRE.sh -h
    ``` 

    Let's set up an example: we want to know the current state-of-the-art regarding chemically-induced colitis in melanoma patients undergoing checkpoint-inhibitors therapy. Our ENQUIRE job might then look something like

    ```bash
    ./code/ENQUIRE.sh -t ICI_and_Colitis -i input/test_input/pmid-ICI_and_Colitis.txt
    ```
    Where all the other parameters described in the `Help` message of `ENQUIRE.sh` are set to default values. The passing of the parameters could be easen by using the `ENQUIRE_config.txt` file that resides in the `/input` subdirectory: the left hand side of each variable assignment must be kept unchanged, while the right hand side can be tweaked according to one's needs. Additional information on the parameters are given in `ENQUIRE_flowchart.png`. Then, the program can be launched by running:

    ```bash
    ./code/ENQUIRE.sh -f input/ENQUIRE_config.txt
    ```
</details>

<details><summary>EXECUTING POST-HOC ANALYSES</summary> 

#### Context-aware gene set annotation 
- Activate `ENQUIRE` environment and run `Rscript code/context_aware_gene_sets.R` from the `ENQUIRE` directory to perform automatic annotation of gene sets, using ENQUIRE-generated, Gene/MeSH edge and node tables and Fuzzy-C-Means (FCM). See the original manuscript for further information.

```
Usage: Rscript code/context_aware_gene_sets.R [options]

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

- You can use the exemplary output files contained in `tmp-FerroDrugTherapy` to test the script:
```bash
Rscript code/context_aware_gene_sets.R -e tmp-FerroDrugTherapy/FerroDrugTherapy/FerroDrugTherapy_Complete_edges_table_subgraph.tsv -n tmp-FerroDrugTherapy/FerroDrugTherapy/FerroDrugTherapy_Complete_nodes_table_subgraph.tsv
```
Please note that the script might last quite long, due to the FCM algorithm.

#### Context-aware pathway enrichment analysis
- Activate `ENQUIRE` environment and run `Rscript code/context_aware_pathway_enrichment.R` from the `ENQUIRE` directory to perform topology-based, pathway enrichment analysis using SANTA and STRING's *H. sapiens*, physical PPI network, using ENQUIRE-generated, gene-gene edge table. See the original manuscript for further information.

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

- You can use the exemplary output files contained in `tmp-FerroDrugTherapy` to test the script:
```bash
Rscript code/context_aware_pathway_enrichment.R -e tmp-FerroDrugTherapy/FerroDrugTherapy/FerroDrugTherapy_Genes_edges_table_subgraph.tsv -c 8 -s 30
```
Please note that the script might last quite long, and it benefits from a high performance computer, if available. 

 
</details>

<details><summary> FREQUENTLY ENCOUNTERED ERRORS </summary>

- It is possible that the R modules fail to proceed because of an error with `pandoc`. If you haven't installed it yet, [proceed to install pandoc](https://pandoc.org/installing.html) and [add pandoc to your PATH environment variable](https://github.com/rstudio/rmarkdown/issues/851).

- Test the command `awk '/MemAvailable/ {print $2}' /proc/meminfo` on your command line: this is the way ENQUIRE checks the available RAM on Linux systems, in order to avoid overflows.  Make sure `awk` is installed on your system. If you witness a non-awk related issue, contact Luca with information on your system and possible solutions to alternatively track the available memory on your OS.

- When computing large networks, an error related to the default `Stack Size` can potentially appear, especially when running R scripts, such as `Error: C stack usage is too close to the limit`. In this case, one shall set a higher stacksize to allow the script to complete, via 

    ```
    ulimit -s N 
    ```    
    Where `N` shall be a size expressed in Kb to set as the maximum stack size. You could first check the number returned by `Cstack_info()` in an active R shell. You can read more about the issue [here](https://stackoverflow.com/questions/14719349/error-c-stack-usage-is-too-close-to-the-limit) and [here](https://rdrr.io/r/base/Cstack_info.html). 

- When running the Fuzzy-C-Means Gene Sets clustering, We observed difficulties and/or failures in installing some of the R packages. In particular, under **Arch Linux** distributions, an *ad hoc* installation of the dependency `V8` might be needed. Please read [here](https://github.com/jeroen/V8/issues/80) about the issue, where it is suggested to install the [Arch Linux specific `v8-r` package](https://aur.archlinux.org/packages/v8-r). This shall allow the script `context_aware_gene_sets.R` to successfully install all R libraries. We recommend trying installing from source with an *ad hoc*, R-version-specific versions of the involved packages, in an R interactive shell.
</details>

<details><summary>IMPORTANT INFORMATION ON PUBMED ACCESSIBILITY</summary>
As of 21.11.22, [important changes](https://www.nlm.nih.gov/pubs/techbull/so22/so22_updated_pubmed_e_utilities.html) have been applied to NCBI's e-utilities. In particular, it is now impossible to stream all records exceeding 10,000 PMIDs from any particular query to the PubMed database. This required to redesign the use of the e-utilities. While it's overall functionality was still preserved, we cannot guarantee the retrieval of all matching records, if the network-based queries obtained by intersecting relevant entities match more than 10,000 records (typically, this is a rare event when intersecting at least 4 distinct entities).
</details>

<details><summary>EXPLANATION OF THE OUTPUT DATA STRUCTURE</summary>

- Provided a recognisable "tag" has been passed to textmining algorithm, a typical output would produce a folder `tmp-tag`, which in turn contains as many subdirectories as the number of steps/iterations performed. For example, if the algorithm constructed 
    
    1. Raw Gene-Mesh network from the original set of papers;
    2. One query expansion while the Gene-Mesh network was not complete yet;
    3. One query-expansion while the Gene-Mesh network was complete, but not the Gene-Gene one;

    Then there will be three subfolders, namely `tag`, `tag_iteration1`, `tag_subgraph_expansion2`. The counter attached to folders and file names records the subsequent attempts to the reconstruction of co-occurence networks. Typically, within each of these sub-folders three pairs of edges and nodes tables can be found corresponding to the respective "Complete" (Gene-Mesh), "Gene" and "Mesh" networks for each iterations (TSV files). These files can be easily imported in Cytoscape or similar graph visualization tools.
    
    Whenever it wasn't possible to obtain one or more of the aforementioned networks, the pipeline should print a message with information on the most meaningful files to look at. It is worth mentioning that the file "tag...Complete_literature_links.tsv" within each subfolder allows fast retrieval of specific edge-associated papers by means of encoded hyperlinks. The batch of queries that were tested in each iteration is stored in "tag...ordered_queries.tsv" within each respective subfolder, with the number of columns corresponding to the `a` attempts at connecting any two communities. Additional meta-data can be explored under the `data/` subfolder. 
    
    Furthemore, under `tmp-tag`, the file `source_pmids.txt` contains all the inspected articles for the given job, which can also be consulted specifically for each iteration under `tmp-tag/efetch_inputs`.
   
    Please contact Luca for any clarification on the purposes of any file.
- NEW: interactive .html networks
    It is now possible to visually inspect Gene-MeSH networks and the reduced networks of entities participating in cliques in two .html files, respectively stored within each iteration's subfolder as "tag...interactive_Gene-MeSH_Network.html" and "tag...interactive_Cliques_Network.html". An exemplary file can be found in the main repository, in the context of a case study of signal transduction pathways in uveal melanoma. 
</details>
