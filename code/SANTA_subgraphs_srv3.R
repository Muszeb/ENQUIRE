#!/home/musellla/miniconda3/envs/R363/lib/R/bin/Rscript
# SANTA RANKINGS BASED ON JACCARD DISTANCE-DERIVED NODE AND EDGE ATTRIBUTES
print("##########################################################################")
print("SANTA: NODE RANKINGS BASED ON JACCARD DISTANCE-DERIVED NODE AND EDGE ATTRIBUTES")
#
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager",repos="https://ftp.fau.de/cran/",lib.loc=libloc)
# CHECK FOR DEPENDENCIES TO BE INSTALLED
# set library location
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("At least two arguments must be supplied (working dir and tag)", call.=FALSE)
} else{
  wd=args[1]
  tag=args[2]
  if (length(args)>=3 & args[3] != "gene" & is.na(as.integer(args[3]))) {
    libloc=args[3] #"/home/musellla/miniconda3/envs/bioinfo/lib/R/library"    
  }
  else if (length(args)>=3 & !(is.na(as.integer(args[3])))) {
    A=as.integer(args[3]) + 1 # used to check for ties within the A+1 highest ranked nodes 
  }
  else if (length(args)>=4 & !(is.na(as.integer(args[4])))) {
    A=as.integer(args[4]) + 1 # used to check for ties within the A+1 highest ranked nodes 
  }
}
if (isFALSE(exists("libloc"))) {
  libloc=.libPaths()[1]  
}
print(c("library in use:",libloc))
#libloc="/home/musellla/miniconda3/envs/bioinfo/lib/R/library"
#
list.of.packages <- c("igraph","BiocManager","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages(lib.loc=libloc)[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc) # FAU mirror
#
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager",repos="https://ftp.fau.de/cran/",lib.loc=libloc)
if(!(c("SANTA") %in% installed.packages(lib.loc=libloc,force=T)[,"Package"])) BiocManager::install("SANTA",lib.loc=libloc,destdir = libloc, force=T) # FAU mirror
#
library(igraph,lib.loc=libloc)
library(utils)
library(SANTA,quietly = T,lib.loc=libloc)
#
print(.libPaths())
#
#### LOAD DATA ####
if ("gene" %in% args){
  #print(paste(wd,"data/gene-subgraphs/",sep = '',collapse = ''))
  setwd(paste(wd,"data/gene-subgraphs/",sep = '',collapse = ''))
} else {
  setwd(paste(wd,"data/subgraphs/",sep = '',collapse = ''))
}
#lib.loc=libloc)
edge_files = unique(grep("edges\\w+",list.files(),value = T)) # list.files() %>% str_subset(pattern = paste(tag,"_edges\\w+",sep = '',collapse = ''))
node_files = unique(grep("nodes\\w+",list.files(),value = T)) #list.files() %>% str_subset(pattern = "nodes\\w+")
edge_files=sort(edge_files);node_files=sort(node_files)
files=list()
for (i in seq(1,length(node_files))){
  files[[node_files[i]]]=c(node_files[i],edge_files[i])  
}
#
# SANTAsub=function(v,args=args){
#   #### LOAD ####
#   nodes=read.delim(v[1],header = TRUE,sep = '\t')#,col.names = c('GENE','MESH','O-R','p2s','pl','pg'))
#   edges=read.delim(v[2],header = TRUE,sep = '\t')#,col.names = c('GENE','MESH','O-R','p2s','pl','pg')) 
#   nodes=nodes[order(nodes$category),]
#   nodes[,-c(1,2)]<-t(apply(nodes[,-c(1,2)],1,signif,digits=3))
#   edges[,c(1,2)]=t(apply(edges[,c(1,2)], 1, sort))
#   edges=edges[order(edges$kw1,edges$kw2),]
#   rownames(edges) <- NULL;rownames(nodes) <- NULL
#   #
#   if (length(colnames(nodes))>6){
#     colnames(nodes)=c('label','category','node_weight','cls_cnt','btw_cnt','strength',colnames(nodes)[seq(7,length(colnames(nodes)))])  
#   } else{
#     colnames(nodes)=c('label','category','node_weight','cls_cnt','btw_cnt','strength')
#   }
#   #,'qscore','Krank','ranking')
#   ##### CREATE GRAPH #####
#   net=graph_from_data_frame(edges, directed = F, vertices = nodes)
#   # assign weights to edges
#   # apply -log10 scale 
#   #
#   net <- set.edge.attribute(net, "weight", index=E(net), edges$inv_edge_weight)#signif(-log10(edges$edge_weight),digits=4)) #signif(1/(2^(edges$R1C1-1)),digits=4)) #signif(-log10(edges$edge_weight),digits=4))  # the smaller, the stronger
#   #### DEFINE SCORE AND APPLY SANTA ####
#   # SUBSTITUTE BTW_CNT WITH ECCENTRICITY 
#   ecc=eccentricity(net)
#   # transform ecc into an increasing function
#   ecc=max(ecc)-ecc
#   net=set.graph.attribute(net,"eccentricity",ecc)
#   nodes$btw_cnt <- ecc[match(nodes$label, names(ecc))]
#   #
#   # maxB=max(nodes$btw_cnt) minB=0
#   # #get ECDF of betweenness centrality, cls cent. and strength
#   distB=ecdf(sort(nodes$btw_cnt)) 
#   distC=ecdf(sort(nodes$cls_cnt))
#   distS=ecdf(sort(nodes$strength))
#   #
#   score=function(l,B,C,S){  #,maxB=maxB,minB=0){#function(l,P,C,B,S){
#     c=(distB(B)*distC(C)*distS(S)) # percentiles product
#     if (c==1){
#       c=0.99 # impede negative or inf values
#     }
#     c=-1/log10(c)
#     #c=distB(B)*distC(C)*distS(S)
#     #c=c(-1/log10((P*C*B*S))) #M=max(B) B=(B-minB)/(maxB-minB) c=S/(2-B)
#     names(c)=l
#     return(c)
#   }
#   qscore <- signif(mapply(score,nodes$label,nodes$btw_cnt,nodes$cls_cnt,nodes$strength),digits = 4)
#   nodes$qscore <- qscore[match(nodes$label, names(qscore))]
#   #
#   net=set_vertex_attr(net, "qscore", index = V(net), nodes$qscore)
#   grank=signif(Knode(net, dist.method=c("shortest.paths"), vertex.attr="qscore", edge.attr="weight", correct.factor = 1, nsteps=10000, B=NULL, verbose=F),digits = 4)
#   # match named rank to a new column in the df
#   nodes$Krank <- grank[match(nodes$label, names(grank))]
#   # include a "true" rank as a scale as well
#   grank=grank[order(grank,names(grank),decreasing = T)]#sort(grank,decreasing = T)
#   if (length(unique(grank[1:A])) < length(grank[1:A])){
#     print(sprintf("### WARNING: there are 2 or more ties within the %i highest-ranked nodes in subgraph %s ###",A,tail(strsplit(v[1],'[_.]')[[1]],n = 2)[1]))
#   }
#   rankings=seq(1,length(grank))
#   names(rankings)=names(grank)
#   nodes$ranking=rankings[match(nodes$label, names(rankings))]
#   return(nodes)
# }
SANTAsub=function(v,args=args){
  #### LOAD ####
  Nodes=read.delim(v[1],header = TRUE,sep = '\t')#,col.names = c('GENE','MESH','O-R','p2s','pl','pg'))
  edges=read.delim(v[2],header = TRUE,sep = '\t')#,col.names = c('GENE','MESH','O-R','p2s','pl','pg')) 
  Nodes=Nodes[order(Nodes$KW),]
  #Nodes[,-c(1,2)]<-t(apply(Nodes[,-c(1,2)],1,signif,digits=3))
  edges[,c(1,2)]=t(apply(edges[,c(1,2)], 1, sort))
  edges=edges[order(edges$kw1,edges$kw2),]
  rownames(edges) <- NULL;rownames(Nodes) <- NULL
  #
  print(head(edges))
  # if (length(colnames(Nodes))>6){
  #   colnames(Nodes)=c('KW','category','node_weight','cls_cnt','btw_cnt','strength',colnames(Nodes)[seq(7,length(colnames(Nodes)))])  
  # } else{
  #   colnames(Nodes)=c('KW','category','node_weight','cls_cnt','btw_cnt','strength')
  # }
  #,'qscore','Krank','ranking')
  ### initialize columns
  Nodes$btw_cnt=NA
  Nodes$cls_cnt=NA
  Nodes$strength=NA
  Nodes$qscore=NA 
  Nodes$Knode=NA
  Nodes$ranking=NA
  ##### CREATE GRAPH #####
  Net=graph_from_data_frame(edges, directed = F, vertices = Nodes)
  Isolated = which(degree(Net)==0)
  Net = delete.vertices(Net, Isolated)
  # assign weights to edges
  # apply -log10 scale 
  #
  edges$inv_edge_weight[edges$inv_edge_weight == 0] <- 1e-9 # vaules of 0 throw an error in Knet
  Net <- set.edge.attribute(Net, "weight", index=E(Net), edges$inv_edge_weight) #signif(1/(2^(edges$R1C1-1)),digits=4)) #signif(-log10(edges$edge_weight),digits=4))  # the smaller, the stronger
  Net <- set.edge.attribute(Net, "R1C1", index=E(Net), edges$R1C1) #signif(1/(2^(edges$R1C1-1)),digits=4)) #signif(-log10(edges$edge_weight),digits=4))  # the smaller, the stronger
  
  #### RUN THE ANALYSIS COMPONENT-WISE ####
  #### DEFINE SCORE AND APPLY SANTA ####
  CC=components(Net)
  for (c in seq(1,CC$no)){
    c=names(CC$membership[CC$membership==c])
    net=induced_subgraph(Net,c)
    nodes=Nodes[Nodes$KW %in% c,]
    # SUBSTITUTE BTW_CNT WITH ECCENTRICITY 
    ecc=eccentricity(net)
    # transform ecc into an increasing function
    ecc=max(ecc)-ecc
    net=set.graph.attribute(net,"eccentricity",ecc)
    nodes$btw_cnt <- ecc[match(nodes$KW, names(ecc))]
    #
    cls=closeness(net,normalized = T)
    nodes$cls_cnt <- cls[match(nodes$KW, names(cls))]
    #
    stre=strength(net,weights = E(net)$R1C1)
    nodes$strength <- stre[match(nodes$KW, names(stre))]
    #
    # get ECDF of betweenness centrality, cls cent. and strength
    distB=ecdf(sort(nodes$btw_cnt)) 
    distC=ecdf(sort(nodes$cls_cnt))
    distS=ecdf(sort(nodes$strength))
    #
    score=function(l,B,C,S){  #,maxB=maxB,minB=0){#function(l,P,C,B,S){
      c=(distB(B)*distC(C)*distS(S)) # percentiles product
      if (c==1){
        c=0.9999 # impede negative or inf values
      }
      c=-1/log10(c)
      #c=distB(B)*distC(C)*distS(S)
      #c=c(-1/log10((P*C*B*S))) #M=max(B) B=(B-minB)/(maxB-minB) c=S/(2-B)
      names(c)=l
      return(c)
    }
    qscore <- signif(mapply(score,nodes$KW,nodes$btw_cnt,nodes$cls_cnt,nodes$strength),digits = 4)
    #nodes$qscore <- qscore[match(nodes$KW, names(qscore))]
    #
    net=set_vertex_attr(net, "qscore", index = V(net), qscore)
    grank=Knode(net, dist.method=c("shortest.paths"), vertex.attr="qscore", edge.attr="weight", correct.factor = 1, nsteps=10000, B=NULL, verbose=F)
    # map metric to complete nodes table
    Nodes[match(names(ecc),Nodes$KW),]$btw_cnt <- ecc
    Nodes[match(names(cls),Nodes$KW),]$cls_cnt <- cls
    Nodes[match(names(stre),Nodes$KW),]$strength <- stre
    Nodes[match(names(qscore),Nodes$KW),]$qscore <- qscore
    Nodes[match(names(grank),Nodes$KW),]$Knode <- grank
    grank=sort(grank,decreasing = T)
    rankings=seq(1,length(grank))
    names(rankings)=names(grank)
    Nodes[match(names(rankings),Nodes$KW),]$ranking <- rankings   
  }
  return(Nodes)
}
t=lapply(files,SANTAsub)

for (i in names(t)){
  #print(i)
  df=t[[i]]
  #print(colnames(df))
  write.table(df,file = i,row.names = F,col.names = T,quote = F,sep = '\t')
  #print(paste(i,'.tsv',sep = '',collapse = ''))
}
print("SANTA: NODE RANKINGS ASSIGNED AND STORED WITHIN THE NODES TABLE")
print("##########################################################################")