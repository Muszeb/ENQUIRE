
library('DBI')
library('dplyr')
library('biomaRt')
library('readr')
library('epanetReader')

mouse_symbols <- read_delim(file='input/mgi_sym_to_aliases.csv')

mouse_genes=mouse_symbols$symbol
mouse_genes=unique(mouse_genes)

### Basic function to convert mouse to human gene names
musGenes <- c("Hmmr")#, "Tlx3", "Cpeb4")
convertMouseGeneList <- function(x){
  x=unique(x)
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = 'mgi_symbol', filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=FALSE)
  #humanx <- unique(genesV2[, 2)
  genesV2 = genesV2[!(duplicated(genesV2[,1])),]
  # Print the first 6 genes found to the screen
  #print(head(humanx))
  return(genesV2)
}

human_mouse=read_delim('input/HOM_MouseHumanSequence.rpt.1')


max <- 10
x <- seq_along(mouse_genes)
split_genes <- split(mouse_genes, ceiling(x/max))

mouse_to_human=lapply(split_genes,convertMouseGeneList)
convertMouseGeneList(musGenes)

# # Basic function to convert human to mouse gene names
# convertHumanGeneList <- function(x){
#   
#   require("biomaRt")
#   human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#   mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#   
#   genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
#   
#   humanx <- unique(genesV2[, 2])
#   
#   # Print the first 6 genes found to the screen
#   print(head(humanx))
#   return(humanx)
# }
# # using retrieved mouse genes as input
# genes <- convertMouseGeneList(musGenes)
# write.table(tibble(genes),'Mouse_TAMs_Textmined_genes_converted_to_Homo.txt', sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
