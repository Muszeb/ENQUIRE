library('DBI')
library('dplyr')
library('biomaRt')
library('readr')
library('epanetReader')

mouse_to_human=read.delim('input/Mouse_to_Human_mart_export.txt')
mouse_to_human=data.frame(mouse_to_human$Human.gene.name,mouse_to_human$Gene.name)
colnames(mouse_to_human)=c("human_genes","mouse_genes")                          
mouse_to_human=mouse_to_human[mouse_to_human$mouse_genes != "",]
mouse_to_human=mouse_to_human[mouse_to_human$human_genes != "",]
mouse_to_human=mouse_to_human[!(duplicated(mouse_to_human)),]

write.table(mouse_to_human,file="input/mouse_to_human_mart_export_cleaned.tsv", sep='\t', col.names = c("human_genes","mouse_genes"), row.names = FALSE, quote = FALSE)
