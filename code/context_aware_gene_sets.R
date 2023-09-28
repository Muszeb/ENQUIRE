#
if (isFALSE(exists("libloc"))) {
  libloc=.libPaths()[1]  
}
#
options(Ncpus = 4L)
print(c("library in use:",libloc))
### install missing libraries ###
print("check for missing, required packages...")
#install.packages(c(),repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc)
###
# factoextra and randomcolor gridBase excluded
if (!("BiocManager" %in% installed.packages(lib.loc=libloc)[,"Package"])){
  install.packages("BiocManager",repos="https://cloud.r-project.org/",lib=libloc,destdir = libloc)
}
list.of.packages <- c("optparse","lsa","RColorBrewer","openxlsx","igraph","BiocManager","stringr","ggplot2","dplyr","dynamicTreeCut","ppclust")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages(lib.loc=libloc)[,"Package"])]
if (length(new.packages)>=1){  
  BiocManager::install(new.packages,lib=libloc,destdir = libloc) # FAU mirror
}
### handle failed installations ###
failed.packages <- new.packages[!(new.packages %in% installed.packages(lib.loc=libloc)[,"Package"])]
if (length(failed.packages)>0){
  print("FAILED PACKAGES:")
  print(failed.packages)
  print("##### if the problem persists, try installing these packages into an interactive R shell with ad hoc package versions #####")
}
# create parser object
library(optparse)
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
option_list <- list( 
    make_option(c("-w", "--directory"), action="store", default=getwd(),
    help="Output directory [default to current working directory]",dest="wd",
    type="character",metavar="path"),
    make_option(c("-e", "--edgetable"),action="store", type="character", default=NULL, 
    help="Path to an ENQUIRE-generated, Gene/MeSH edge table file (required)",
    metavar="path",dest="test_net"),
    make_option(c("-n", "--nodetable"),action="store", type="character", default=NULL, 
    help="Path to an ENQUIRE-generated, Gene/MeSH node table file (required)",
    metavar="path",dest="test_nodes"),
    make_option(c("-t", "--tag"),action="store", type="character", default="ENQUIRE", 
    help="tag prefix (default to 'ENQUIRE')",
    metavar="tag",dest="tag"),
    make_option(c("-d", "--membdeg"),action="store", type="numeric", default=.05, 
    help="minimal membership degree for gene-to-cluster association (default: 0.05), range [0-1]",
    metavar="parameter",dest="d"),
    make_option(c("-s", "--setsize"),action="store", type="numeric", default=2, 
    help="minimal gene set size (default: 2)",
    metavar="parameter",dest="s"))

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))
if (length(opt)<5){
  print("ERROR: missing inputs!")
  print("Do 'Rscript context_aware_gene_sets.R -h' for printing help")
  stop()
}
print("parsed options:")
print(str(opt))
print("#####")
###
library(openxlsx)
library(ppclust)
library(igraph)
library(dplyr)
library(ggplot2)
library(dynamicTreeCut)
library(stringr)
library(RColorBrewer)
#### 
test_net <- read.delim(opt$test_net,stringsAsFactors = F)
test_nodes <- read.delim(opt$test_nodes,stringsAsFactors = F)
#
net=graph_from_data_frame(test_net,directed = F,vertices = test_nodes)
net=delete.vertices(net,which(degree(net)==0))
##
## use invlogweighted similarity ####
##
sims=matrix(similarity(net,method='invlogweighted'),
            ncol=vcount(net),
            nrow=vcount(net),
            dimnames = list('colnames'=V(net)$name,'rownames'=V(net)$name))

###
### use ward clustering to simplify genes sets - at least 2 entities per set
###
### scale to reduce the impact of variance on prioritizing directions
x=dist(scale(sims),method = 'euclidean')
###
wardtree=hclust(x,method = 'ward.D2')
#wardtree$height <- round(wardtree$height, 6) 
clustmesh=setNames(dynamicTreeCut::cutreeDynamic(wardtree,minClusterSize = 2,distM = as.matrix(x)),wardtree$labels)
print("WARD AND DYNAMIC TREE CUT RESULTS:")
print(table(clustmesh))
### extract condensed centers by specifying memberships from ward/treecut ###
### also known as sharding ###
shards=tapply(names(clustmesh),factor(unname(clustmesh)),function(rows)(apply(sims[rows,],2,mean)))
shards=do.call('rbind',shards)
### use ward for initial memberships ###
uStart=tapply(names(clustmesh)[order(unname(clustmesh))],factor(sort(unname(clustmesh))),function(rows){
  return(matrix(0.0,nrow = length(rows),ncol=length(levels(factor(unname(clustmesh)))),
                dimnames = list(rows,levels(factor(unname(clustmesh))))))})
#
for (i in names(uStart)){
  uStart[[i]][,i]=rep(1,nrow(uStart[[i]]))
}
uStart=do.call('rbind',uStart)
### FCM ###
x=as.matrix(x)
print("perform FCM with Ward-derived initial centers...")
res.fcm <- ppclust::fcm(x,stand = F,memberships = uStart,centers = shards,m=2,fixcent = F,nstart=3,iter.max = 1E6,con.val=1e-9,numseed = 2202)
#
genes=test_nodes[test_nodes[,2,drop=T]=='GENE',1,drop=T]
print(genes)
#
print("...done, produce output data...")
#### use PCA then custom plot ####
pcadat=prcomp(x,center = T)
###
### use defuzzyfied clusters to draw ellipses 
#res.fcm$cluster=setNames(res.fcm$cluster,rownames(sims))
clusts=setNames(res.fcm$cluster,rownames(sims))#setNames(unname(res.fcm$cluster),rownames(sims)[match(names(res.fcm$cluster),rownames(sims))])
pcashard=tapply(names(clusts),factor(unname(clusts)),function(v)(as.data.frame(matrix(pcadat$x[v,c('PC1','PC2')],
                                                                                      ncol = 2,
                                                                                      nrow = length(v),
                                                                                      dimnames = list(rownames=v,colnames=c('PC1','PC2'))))))
pcashardf=bind_rows(pcashard)
pcashardf$Cluster=clusts[match(rownames(pcashardf),names(clusts))] #unlist(sapply(names(pcashard),function(n)(rep(n,nrow(pcashard[[n]])))))
pcashard=pcashardf
#pcashard=cbind(data.frame('Cluster'=unlist(sapply(names(pcashard),function(n)(rep(n,nrow(pcashard[[n]]))))),row.names = names(pcashard)),
#               bind_rows(pcashard))
### gene set parameters ###
d=opt$d
s=opt$s
### genes set of at least 2 - given minimum memberships to cluster ###
okclust=as.character(seq(1,ncol(res.fcm$u))[apply(res.fcm$u[genes,],2,function(v)(sum(as.numeric(v>=d))>=s))])
###
pcashard$Cluster=factor(as.numeric(pcashard$Cluster),levels = sort(as.numeric(unique(pcashard$Cluster))))
##### implement cluster labels ####
#pcashard$sorting=test_nodes[match(rownames(pcashard),test_nodes$KW),ncol(test_nodes)-2,drop=T]
ttt=function(row,clus)(res.fcm$u[row,clus])
pcashard$MembDeg=mapply(FUN = function(row,clus)(res.fcm$u[row,clus]),
                        rownames(pcashard),as.numeric(as.character(pcashard$Cluster)))
#pcashard$ranklab=mapply(FUN = function(q,u)(q*u),pcashard$sorting,pcashard$membdeg)
#
# top meshs by closeness to centroid define descriptions #
vectbyclus=tapply(as.list(as.data.frame(t(x[!(rownames(x) %in% genes),]))),
                  factor(res.fcm$cluster[!(rownames(x) %in% genes)]),function(l)(l))# as.list(as.data.frame(t(res.fcm$v)))
#
vectmin=lapply(as.list(as.data.frame(t(res.fcm$v))),function(clus){
  eucres=apply(x[!(rownames(x) %in% genes),],1,function(v)(sqrt(sum((clus - v)^2))))
  return(names(eucres[order(eucres)[1:4]]))})
#
names(vectmin)=as.numeric(str_extract(names(vectmin),pattern = '[0-9]+'))
#
labs=unlist(lapply(vectmin,function(x)(paste(str_trunc(str_remove_all(x,"\\/.+$"),width = 18), sep = '',collapse=', '))))
labslong=unlist(lapply(vectmin,function(x)(paste(x,sep = '',collapse=', '))))
###
pcashard$lab=labs[match(pcashard$Cluster,names(labs))]
pcashard$lablong=labslong[match(pcashard$Cluster,names(labslong))]
#### save also complete mesh-cluster associations ####
fullmesh2cl=data.frame(cl=unname(res.fcm$cluster),en=rownames(res.fcm$u))
#fullmesh2cl=fullmesh2cl[!(fullmesh2cl$en %in% rownames(testmat)), ]
fullmeshs=tapply(fullmesh2cl$en,factor(fullmesh2cl$cl),function(v)(paste(v,sep=', ',collapse = ' | ')))
pcashard$Category=test_nodes[match(rownames(pcashard),test_nodes$KW),2,drop=T]
clust2name=setNames(pcashard$Cluster,pcashard$lab)
clust2name=clust2name[!(duplicated(names(clust2name)))]
####
### PIE CHART REPRESENTATION ###
props=lapply(rownames(res.fcm$u),function(r)(setNames(res.fcm$u[r,,drop=T],colnames(res.fcm$u))))
names(props)=rownames(res.fcm$u)

### account for vectors of just 0s
if (all(sapply(props,function(v)(all(v<=0.05))))){
  print("WARNING: ALL PREDICTED MEMBERSHIP DEGREES ARE < 5%")
  print("The clusters are highly fuzzy")
  print("proportions will be kept intact")
}else{
  print('simplify memberships based on a minimum of 5% memberships')
  props=lapply(props,function(v)(ifelse(v >= d,v,0)))
  props[sapply(props,function(v)(all(v < d)))]=lapply(props[sapply(props,function(v)(all(v < d)))],function(v)(rep(1/length(v),length(v))))  
}
#
#### GENERATE OUTPUT TABLES ####
library(openxlsx)
wb <- createWorkbook()
# 1)= genes x cluster name membership degree #
addWorksheet(wb, "GenesXClusterMembDeg",gridLines = T)
gcmat=res.fcm$u[match(genes,rownames(res.fcm$u)),]
writeDataTable(wb, "GenesXClusterMembDeg", as.data.frame(gcmat), rowNames = T)
# 2) gene sets according to k/d
addWorksheet(wb, "FuzzyGeneSets",gridLines = T)
gcl=apply(gcmat,2,function(cls)(rownames(gcmat)[cls>=d]))
gcl=gcl[which(lengths(gcl)>=s)]
ml=max(lengths(gcl))
gcl=lapply(gcl,function(v)(c(v,rep('',ml-length(v)))))
gcl=rbind(labslong[as.numeric(str_match(names(gcl),"[0-9]+")[,1])],as.data.frame(gcl))
rownames(gcl)=c('MeSH Centroids (4)',seq(1,nrow(gcl)-1))
writeDataTable(wb, "FuzzyGeneSets", gcl, rowNames = T)
# 3) crisp clusters 
addWorksheet(wb, "PCA-CrispClusters",gridLines = T)
writeDataTable(wb, "PCA-CrispClusters", pcashard, rowNames = T)
# save!
setwd(opt$wd)
xlsxfile = sprintf("%s_context_aware_gene_sets.xlsx",opt$tag)
saveWorkbook(wb, xlsxfile, overwrite = T) 