options(Ncpus = 8L)
###
if (isFALSE(exists("libloc"))) {
  libloc=.libPaths()[1]  
}
#
print(c("library in use:",libloc))
### install missing libraries ###
print("check for missing, required packages...")
list.of.packages <- c("optparse","lsa","openxlsx","RColorBrewer","dynamicTreeCut","ppclust","RColorBrewer","qgraph","dynamicTreeCut","fclust","svglite","gridBase","igraph","BiocManager","tidyr","parallel","reshape2","stringr","ggplot2","snow","dplyr","networkD3","htmlwidgets","data.table","poolr","VarianceGamma","magrittr","pkgconfig")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages(lib.loc=libloc)[,"Package"])]
if(length(new.packages)){
  options(Ncpus = 8)
  install.packages(new.packages,repos="https://cloud.r-project.org/",lib=libloc,destdir = libloc) # FAU mirror
}
#if(!(c("SANTA") %in% installed.packages(lib.loc=libloc,force=T)[,"Package"])) BiocManager::install("SANTA",lib.loc=libloc,destdir = libloc, force=T) # FAU mirror
# if(!(c("SANTA") %in% installed.packages(lib.loc=libloc,force=T)[,"Package"])){
#   options(Ncpus = 8)
if(!(c("SANTA") %in% installed.packages(lib.loc=libloc,force=T)[,"Package"])) BiocManager::install("SANTA",lib.loc=libloc,destdir = libloc, force=T) # FAU mirror
###