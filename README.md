# ENQUIRE
<ins>E</ins>xpanding <ins>N</ins>etworks by <ins>Q</ins>uerying <ins>U</ins>nexpectedly <ins>I</ins>nter-<ins>R</ins>elated <ins>E</ins>ntities

## ABSTRACT 
The accelerating growth in scientific literature is overwhelming our capacity to manually distil complex phenomena like molecular networks linked to diseases. Moreover, confirmation biases in search engines and databases influence the interpretation of facts and the generation of hypotheses. ENQUIRE (Expanding Networks by Querying Unexpectedly Inter-Related Entities) offers an alternative to manual literature curation and database mining to study complex biomedical phenomena. ENQUIRE generates a co-occurrence network of genes and biomedical ontologies (MeSH) using a corpus of publications as input. The algorithm iteratively reconstructs and expands the network by generating PubMed queries from significant interrelations, until it obtains an interconnected network of genes and MeSH relevant to the input corpus. We systematically evaluated ENQUIRE’s versatility and found that it generates co-occurrence networks similar to those based on co-expression data and manually annotated databases with a fraction of the time and human resources. Using case studies spanning cancer, cell differentiation and immunity, ENQUIRE proved to identify interlinked genes and enriched pathways unique to each topic, thereby preserving their underlying diversity. ENQUIRE supports biomedical researchers by easing literature annotation, boosting hypothesis formulation and facilitating the identification of molecular targets for subsequent experimentation.

## EXEMPLARY RUN

<details><summary>EXPLANATION OF THE OUTPUT DATA STRUCTURE</summary>

- Provided a recognisable "tag" has been passed to ENQUIRE, a typical output would produce a folder `tmp-tag`, which in turn contains as many subdirectories as the number of steps/iterations performed. For example, if the algorithm performed one initial network reconstruction and two successful literature query steps that yielded that many expanded networks, then there would be three subfolders, namely `tag`, `tag_subgraph_expansion1`, and `tag_subgraph_expansion2`. The counter attached to folders and file names records the subsequent attempts at expanding the initial co-occurence network. Typically, within each of these sub-folders, three pairs of edges and nodes tables can be found corresponding to the respective "Complete" (Gene/Mesh), "Gene" and "Mesh" networks for each iterations (TSV files). These files can be easily imported in Cytoscape or similar graph visualization tools.

-  It is also possible to visually inspect Gene/MeSH networks and the reduced networks containing only cliques by means of two `.html` files, respectively stored within each iteration's subfolder as `tag...interactive_Gene-MeSH_Network.html` and `tag...interactive_Cliques_Network.html`.
      
- Whenever it wasn't possible to obtain one or more of the aforementioned networks, the pipeline should print a message with information on the most meaningful files to look at. It is worth mentioning that the file `tag...Complete_literature_links.tsv` within each subfolder allows for fast retrieval of specific edge-associated papers by means of encoded hyperlinks, as well as subsetting the results by input or ENQUIRE-generated PubMed query. Additional meta-data can be explored under the `data/` subfolder.

- Furthemore, under `tmp-tag`, the file `source_pmids.txt` contains all the inspected articles for the given job, which can also be consulted specifically for each iteration under `tmp-tag/efetch_inputs`. A record of all retrieved PMIDs and query attempts can be found under `tmp-tag/efetch_inputs/`

## CODE 
A full release of ENQUIRE standalone will soon be pushed to this repository, so keep an eye on it! :) 



