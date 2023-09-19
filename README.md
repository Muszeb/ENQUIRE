# ENQUIRE

<details><summary>INITIALIZATION</summary> 

ENQUIRE can currently be installed and run on UNIX systems such as MacOS, Linux and UNIX emulators/virtual machines for Windows such as WSL, CygWin or MSYS2 (installation of UNIX emulators is not included in this guide).

Start by downloading the repository folder stored under the branch "ENQUIRE".

The script requires Python 3.8 or later versions (tested up to 3.9.X) on Linux and UNIX emulators, while Python 3.7 is necessary to run ENQUIRE on MacOS (ad hoc build is under construction).
To successfully run ENQUIRE, we need to successfully install **Python and its required libraries**, **EDirect** and an **R version 3.X**. This section is therefore divided into these free subcapters, as well as an additional **frequently encountered errors** section.  

<details><summary>INSTALL PYTHON AND REQUIRED LIBRARIES </summary> 

<details><summary>THROUGH CONDA (RECOMMENDED)</summary> 

The easiest way to install all required packages is to create a new `conda` environment: the algorithm will therefore only work if such environment is "activated". 
To create a dedicated `conda` environment, first install an Anaconda or Miniconda distribution (e.g.: [from here](https://repo.continuum.io/miniconda/)).
For example, on Linux you could run: 

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

and on MacOS:

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

- Once `conda` has being installed, use the file called `ENQUIRE.yml` that resides in the `input` subdirectory to generate a new environment called "ENQUIRE" (or any other name you like), by running from the main directory of this repository (Windows user may want to first [locate their local disk from a UNIX virtual machine](https://askubuntu.com/questions/943006/how-to-navigate-to-c-drive-in-bash-on-wsl-ubuntu)): 

```
conda env create --name ENQUIRE -f input/ENQUIRE.yml
```

The environment can then be activated by running:

```
conda activate ENQUIRE
```
The environment will have Python 3.8.X and R 3.6.3 installed by default.

- Sometimes, an installation conflict can happen, such that conda does not install some of the required packages. We recommend   running a trial of the installation while redirecting the stdout 2 to a log file (add `2> logfile.log` to the installation command). In case the installation of a package fails, we recommend running an ad-hoc pip installation like `pip install package==x.x.x` within the ENQUIRE environment, choosing a package version as close as possible to those reported in input/ENQUIRE_pip_requirements.txt.
For example, the package `en-ner-jnlpba-md` (check [scispacy's repository](https://allenai.github.io/scispacy/) cannot be installed through `conda`. Then, run the following command to install the last missing packages:   
```
conda activate ENQUIRE
pip install -r input/PackagesNotFound.txt
```
</details>

<details><summary>WITHOUT ANY ENVIRONMENT</summary>

Alternatively, to install all the required Python libraries, install `python >= 3.8` and `pip >= 20.2`, then from the main directory run

```
pip install -r input/ENQUIRE_pip_requirements.txt
```
</details>
</details>

<details><summary> INSTALLING R 3.X </summary>

If You followed the installation through `conda`, R version 3.6.3 should be installed under the `ENQUIRE` environment with its corresponding `Rscript` interpreter. The latter is required to run node ranking functions from the [SANTA](http://bioconductor.org/packages/release/bioc/html/SANTA.html) package. If needed, a second yml file named `R363_env.yml` can be found in the `input/` folder of the main repository. Then, to install `R 3.6.3` separately, simply run 

```
conda env create --name R363 -f input/R363_env.yml
```

Now it is necessary to locate the `Rscript` interpreter. Its file path needs to be passed to ENQUIRE whenever a job is launched. To get the `Rscript` file path, you could run: 
        
```
conda activate ENQUIRE
which Rscript    
```
The printed output should look something like `/home/usr/miniconda3/envs/ENQUIRE/bin/Rscript`. One way to let ENQUIRE know the `RScript` path would be to `export` the file path via (the variable name is important!):

```
export rscript=$(which Rscript)
```

which will let ENQUIRE know the `Rscript` location within the current terminal session. To have this set for any terminal session, write the `rscript` variable into a terminal profile file, such as the `.bashrc`. 
If you don't want to generate the `rscript` variable in your interactive bash shell, just pass `$(which Rscript)` to the `-w` argument of `ENQUIRE.sh`.
Alternatively, copy paste the exact file path and add the line `rscript=/home/usr/miniconda3/envs/R363/bin/Rscript` to your `textmining_config.txt` file (more on that under "Launching the algorithm").


</details>

<details><summary>INSTALL EDIRECT FROM NCBI (E-UTILITIES) </summary> 

As stated in [this document from NCBI](https://www.ncbi.nlm.nih.gov/books/NBK179288/) (Created: April 23, 2013; Last    Update: January 4, 2023.):

EDirect will run on Unix and Macintosh computers, and under the Cygwin or [WSL (Windows 10+ users)](https://ubuntu.com/wsl) Unix-emulation environments on Windows PCs. To install the EDirect software, open a terminal window and execute one of the following two commands:

```
sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
```

```
sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
```

This will download a number of scripts and several precompiled programs into an "edirect" folder in the user's home directory. It may then print an additional command for updating the PATH environment variable in the user's configuration file. The editing instructions will look something like:

```
echo "export PATH=\$PATH:\$HOME/edirect" >> $HOME/.bash_profile
```

As a convenience, the installation process ends by offering to run the PATH update command for you. Answer "y" and press the Return key if you want it run. If the PATH is already set correctly, or if you prefer to make any editing changes manually, just press Return.

Once installation is complete, run:

```
export PATH=${PATH}:${HOME}/edirect
```

to set the PATH for the current terminal session.

</details>

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
</details>

<details><summary>HOW TO EXTRACT LITERATURE FROM REVIEW PAPERS</summary>

- It is now possible to extract the PubMed identifiers of all papers cited in a reading of interest (e.g. a review paper of a particular topic), by means of the python script `efetch_references.py`. From the main repository, type on the command line 

```
python code/efetch_references.py tag ref1 ref2 ref3 ...

```

where `tag` is the name of the plain text output file, while `ref1 ref2 ref3 ...` are the PMIDs of the papers you want to extract the references from. The output will look like the example from the previous section and is therefore ready to be used in the textmining pipeline. 
DISCLAIMER: if the references are not annotated into the Pubmed's API, an error such as 

`File "code/efetch_references.py", line 28, in <module>
    refs+=refparse(p)
  File "code/efetch_references.py", line 20, in refparse
    refs=dpath.get(data,"**/Link") # list of {Id:value} dicts
  File "/home/musellla/miniconda3/envs/wokenv/lib/python3.8/site-packages/dpath/util.py", line 178, in get
    raise KeyError(glob)
KeyError: '**/Link'`

might occur. As a rule of thumb, check the "page navigation" menu on the Pubmed page of the article of interest ([example of a review with non indexed references](https://pubmed.ncbi.nlm.nih.gov/33086849/)). 

</details>

<details><summary>LAUNCHING THE ALGORITHM</summary> 

- After the download, you should see a folder called `tam_textmining`: this is the main directory from which the program is supposed to be run.

- Before running an actual task, it is recommended to inspect the flowchart stored as a PDF in the main repository and the Help section of the code by running (from the main repo directory):
 
    ```
    ./code/ENQUIRE.sh -h
    ``` 

    Let's set up an example: we want to know the current state-of-the-art regarding chemically-induced colitis in melanoma patients undergoing checkpoint-inhibitors therapy. A typical job might look something like

    ```
    ./code/ENQUIRE.sh -p $(pwd)/ -t ICI_and_Colitis -j 6 -i test_input/pmid-ICI_and_Colitis.txt -r 1 -c 4 -a 3 -k 3 -w $rscript
    ```
    but the passing of the parameters could be easen by using the `textmining_config.txt` file that resides in the `/input` subdirectory: the left hand side of each variable assignment must be left unchanged, while the right hand side can be tweaked according to one's needs. Then the program could be launched by running:

    ```
    ./code/ENQUIRE.sh -f input/textmining_config.txt
    ```
At the first execution of the above code using a fresh R distibution, ENQUIRE will launch the installation of all the necessary libraries, which will take several minutes. If the `.libPaths()` object behaves as expected, the process should go smoothly and happen only once.  

</details>

<details><summary> FREQUENTLY ENCOUNTERED ERRORS </summary>

- It is possible that the R modules fail to proceed because of an error with `pandoc`. If you haven't installed it yet, [proceed to install pandoc](https://pandoc.org/installing.html) and [add pandoc to your PATH environment variable](https://github.com/rstudio/rmarkdown/issues/851).

- Test the command `awk '/MemAvailable/ {print $2}' /proc/meminfo` on your command line: this is the way ENQUIRE checks the available RAM on Linux systems, in order to avoid overflows.  Make sure `awk` is installed on your system. If you witness a non-awk related issue, contact Luca with information on your system and possible solutions to alternatively track the available memory on your OS.

- When computing large networks, an error related to the default `Stack Size` can potentially appear, especially when running R scripts, such as `Error: C stack usage is too close to the limit`. In this case, one shall set a higher stacksize to allow the script to complete, via 

    ```
    ulimit -s N 
    ```    
    Where `N` shall be a size expressed in Kb to set as the maximum stack size. You could first check the number returned by `Cstack_info()` in an active R shell. You can read more about the issue [here](https://stackoverflow.com/questions/14719349/error-c-stack-usage-is-too-close-to-the-limit) and [here](https://rdrr.io/r/base/Cstack_info.html). 

- When running the Fuzzy-C-Means Gene Sets clustering, We observed difficulties and/or failures in installing the R packages `factoextra`, `qgraph`, `car` and `randomcoloR` under R 3.6.3. In particular, under **Arch Linux** distributions,`randomcoloR` requires an *ad hoc* installation of the dependency `V8`. Please read [here](https://github.com/jeroen/V8/issues/80) about the issue, where it is suggested to install the [Arch Linux specific `v8-r` package](https://aur.archlinux.org/packages/v8-r). This shall allow the script `fuzzy_clustering_genesets_executive.R` to successfully install all R libraries. For packages `factoextra`, `qgraph` and `car`, an *ad hoc* "installation recipe" that compiles from specific versions of the repositories is executed within the code, hence these packages should hopefully not give any problem.

</details>

<details><summary>NEW: IMPORTANT UPDATES ON PUBMED'S ACCESSIBILITY</summary>
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
