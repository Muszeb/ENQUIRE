#!/home/musellla/miniconda3/envs/R363/lib/R/bin/Rscript
# SANTA RANKINGS BASED ON JACCARD DISTANCE-DERIVED NODE AND EDGE ATTRIBUTES
options(Ncpus = 8L)
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("At least three arguments must be supplied (working dir, tag and ncores)", call.=FALSE)
}else{
  #print(args)
  wd=args[1]
  tag=args[2]
  numargs=as.integer(args[!(is.na(as.integer(args)))])
  #print(numargs)
  if (length(numargs)>=3){
    max_attempts=as.integer(numargs[2])
    combs=as.integer(numargs[3])
    ncores=floor(numargs[1])
    #print(combs)
    A=max_attempts + 1 # used to check for ties within the A+1 highest ranked nodes
  }
  else if (length(numargs)>=2){
    ncores=floor(numargs[1])
    max_attempts=as.integer(numargs[2])
    combs=4
    #print(combs)
    A=max_attempts + 1 # used to check for ties within the A+1 highest ranked nodes
  }
  else{
    ncores=floor(numargs[1])
    max_attempts=3
    combs=4
  }
  # if (length(args)>=3 & args[3] != "incomplete" & is.na(as.integer(args[3]))){
  #   libloc=args[3] #"/home/musellla/miniconda3/envs/bioinfo/lib/R/library"    
  # }
  # else if (length(args)>=3 & !(is.na(as.integer(args[3])))) {
  #   max_attempts=as.integer(args[3])
  #   A=max_attempts + 1 # used to check for ties within the A+1 highest ranked nodes 
  # }
  # else if (length(args)>=4 & !(is.na(as.integer(args[4])))) {
  #   max_attempts=as.integer(args[4])
  #   A=max_attempts + 1 # used to check for ties within the A+1 highest ranked nodes 
  # }
  # else if (length(args)>=5 & !(is.na(as.integer(args[4])))) {
  #   print("it happened")
  #   max_attempts=as.integer(args[4])
  #   combs=as.integer(args[5])
  #   print(combs)
  #   A=max_attempts + 1 # used to check for ties within the A+1 highest ranked nodes 
  # }
}
#
if (isFALSE(exists("libloc"))) {
  libloc=.libPaths()[1]  
}
#
print(c("library in use:",libloc))
### install missing libraries ###
print("check for missing, required packages...")
list.of.packages <- c("igraph","RColorBrewer","BiocManager","tidyr","parallel","reshape2","stringr","ggplot2","snow","dplyr","networkD3","htmlwidgets","data.table","poolr","VarianceGamma","magrittr","pkgconfig")
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
library(utils)
library(SANTA,quietly = T,lib.loc=libloc)
library(stringr,quietly = T,lib.loc=libloc)
library(ggplot2,quietly = T,lib.loc=libloc)
library(parallel,quietly = T,lib.loc=libloc)
library(reshape2,quietly = T,lib.loc=libloc)
library(dplyr)
library(networkD3)
library(htmlwidgets)
library(dplyr)
library(poolr)
library(VarianceGamma)
library(data.table)
library(RColorBrewer)
###
#### LOAD DATA ####
if ("incomplete" %in% args){
  #print(paste(wd,"data/gene-subgraphs/",sep = '',collapse = ''))
  setwd(paste(wd,"data/",sep = '',collapse = ''))
  edge_files = unique(grep("allxall_sig_combinations_bh_to_cytoscape.tsv",list.files(),value = T)) # list.files() %>% str_subset(pattern = paste(tag,"_edges\\w+",sep = '',collapse = ''))
  node_files = unique(grep("_allxall_nodes_table_bh_to_cytoscape.tsv",list.files(),value = T)) #list.files() %>% str_subset(pattern = "nodes\\w+")
} else {
  setwd(wd)
  edge_files = unique(grep("_edges_table",list.files(),value = T)) # list.files() %>% str_subset(pattern = paste(tag,"_edges\\w+",sep = '',collapse = ''))
  node_files = unique(grep("_nodes_table",list.files(),value = T)) #list.files() %>% str_subset(pattern = "nodes\\w+")
}
#
print(node_files)
print(edge_files)
#
edge_files=sort(edge_files);node_files=sort(node_files)
files=list()
for (i in seq(1,length(node_files))){
  files[[node_files[i]]]=c(node_files[i],edge_files[i])  
}
#### KNODE ####
print("##########################################################################")
print("SANTA: NODE RANKINGS BASED ON PRODUCT OF DIRECTIONAL, CONDITIONAL PROBABILITIES")
#
SANTAsub=function(v,args=args){
  #### LOAD ####
  Nodes=read.delim(v[1],header = TRUE,sep = '\t')#,col.names = c('GENE','MESH','O-R','p2s','pl','pg'))
  edges=read.delim(v[2],header = TRUE,sep = '\t')#,col.names = c('GENE','MESH','O-R','p2s','pl','pg')) 
  Nodes=Nodes[order(Nodes$KW),]
  # Nodes[,-c(1,2)]<-t(apply(Nodes[,-c(1,2)],1,signif,digits=3))
  edges[,c(1,2)]=t(apply(edges[,c(1,2)], 1, sort))
  edges=edges[order(edges$kw1,edges$kw2),]
  rownames(edges) <- NULL;rownames(Nodes) <- NULL
  # #
  # if (length(colnames(Nodes))>6){
  #   colnames(Nodes)=c('KW','category','node_weight','cls_cnt','btw_cnt','strength',colnames(Nodes)[seq(7,length(colnames(Nodes)))])  
  # } else{
  #   colnames(Nodes)=c('KW','category','node_weight','cls_cnt','btw_cnt','strength')
  # }
  colnames(Nodes)=c('KW','category','node_weight',colnames(Nodes)[seq(4,length(colnames(Nodes)))])[1:ncol(Nodes)]
  # #,'qscore','Krank','ranking')
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
  edges$inv_edge_weight[edges$inv_edge_weight == 0] <- 1e-15 # vaules of 0 throw an error in Knet
  Net <- set.edge.attribute(Net, "weight", index=E(Net),as.numeric(format(edges$inv_edge_weight,digits=15))) #signif(1/(2^(edges$R1C1-1)),digits=4)) #signif(-log10(edges$edge_weight),digits=4))  # the smaller, the stronger
  Net <- set.edge.attribute(Net, "R1C1", index=E(Net), edges$R1C1) #signif(1/(2^(edges$R1C1-1)),digits=4)) #signif(-log10(edges$edge_weight),digits=4))  # the smaller, the stronger  
  #### RUN THE ANALYSIS COMPONENT-WISE ####
  #### DEFINE SCORE AND APPLY SANTA ####
  CC=components(Net)
  for (c in seq(1,CC$no)){
    print(sprintf("\npercentiles of edge weight, component %i:",as.integer(c)))
    c=names(CC$membership[CC$membership==c])
    net=induced_subgraph(Net,c)
    nodes=Nodes[Nodes$KW %in% c,]
    # SUBSTITUTE BTW_CNT WITH ECCENTRICITY 
    #ecc=eccentricity(net)
    # transform ecc into an increasing function
    #ecc=max(ecc)-ecc
    #net=set.graph.attribute(net,"eccentricity",ecc)    
    print(quantile(E(net)$weight))
    btw=betweenness(net,weights = E(net)$weight)
    nodes$btw_cnt <- btw[match(nodes$KW, names(btw))]  #ecc[match(nodes$KW, names(ecc))]
    #
    if (all(distances(net,V(net),V(net))==0)){
        nodes$cls_cnt <- 0
        cls=nodes$cls_cnt
      }else{
        cls=closeness(net,normalized = T)
        nodes$cls_cnt <- cls[match(nodes$KW, names(cls))] 
      }
    #
    stre=strength(net,weights = as.numeric(format(1-E(net)$weight,digits=15)))
    nodes$strength <- stre[match(nodes$KW, names(stre))]
    #
    # get ECDF of betweenness centrality, cls cent. and strength
    distB=ecdf(sort(nodes$btw_cnt)) 
    #distC=ecdf(sort(nodes$cls_cnt))
    distS=ecdf(sort(nodes$strength))
    #
    score=function(B,S){  #,maxB=maxB,minB=0){#function(l,P,C,B,S){
      c=(distB(B)*distS(S))# # percentiles product
      if (c==0){
        c=1e-15 # impede negative or inf values
      }
      #c=-log10(c)
      #c=distB(B)*distC(C)*distS(S)
      #c=c(-1/log10((P*C*B*S))) #M=max(B) B=(B-minB)/(maxB-minB) c=S/(2-B)
      #names(c)=l
      return(c)
    }
    qscore <- setNames(mapply(score,nodes$btw_cnt,nodes$strength),nodes$KW)
    #nodes$qscore <- qscore[match(nodes$KW, names(qscore))]
    #names(qscore)=stringr::str_remove(names(qscore),"\\.([^.]+)$")
    #print(qscore)
    #nodes$qscore <- qscore[match(nodes$KW, names(qscore))]
    #
    #print(setdiff(names(qscore),Nodes$KW))
    #print(setdiff(Nodes$KW,names(qscore)))
    #
    net=set_vertex_attr(net, "qscore", index = V(net), qscore)
    grank=Knode(net, dist.method=c("shortest.paths"), vertex.attr="qscore", edge.attr="weight", correct.factor = 1, nsteps=10000, B=NULL, verbose=F)
    # map metric to complete nodes table
    Nodes[match(names(btw),Nodes$KW),]$btw_cnt <- btw
    Nodes[match(names(cls),Nodes$KW),]$cls_cnt <- cls
    Nodes[match(names(stre),Nodes$KW),]$strength <- stre
    Nodes[match(names(qscore),Nodes$KW),]$qscore <- qscore
    Nodes[match(names(grank),Nodes$KW),]$Knode <- grank
    grank=sort(grank,decreasing = T)
    rankings=seq(1,length(grank))
    names(rankings)=names(grank)
    Nodes[match(names(rankings),Nodes$KW),]$ranking <- rankings   
  }
  return(list('nodes'=Nodes,'edges'=edges))
}

t=lapply(files,SANTAsub)

for (i in names(t)){
  #print(i)
  df=t[[i]]$nodes
  #print(colnames(df))
  write.table(df,file = i,row.names = F,col.names = T,quote = F,sep = '\t')
  #print(paste(i,'.tsv',sep = '',collapse = ''))
}
print("SANTA: NODE RANKINGS ASSIGNED AND STORED WITHIN THE NODES TABLE")
print("##########################################################################")
#####
print("##########################################################################")
print("SANTA: KNET STATISTICS ON MAXIMAL CLIQUES")
# set working directory to extracted folder
setwd(wd)
# set weight
edges=t[[1]]$edges
edges$weight=edges$inv_edge_weight
nodes=t[[1]]$nodes
#
### GENERATE HTML REPRESANTION IF ALL COLUMNS ARE AVAILABLE ###
#
widnet=function(links,nodes,source,target,nodeid,edge_weight,col_weight,group,nodesize,click=T){
  links$V1=apply(as.matrix(links[,c(source)]),1,function(v)(
    as.numeric(rownames(nodes[nodes[,c(nodeid),drop=T]==v,]))-1))
  links$V2=apply(as.matrix(links[,c(target)]),1,function(v)(
    as.numeric(rownames(nodes[nodes[,c(nodeid),drop=T]==v,]))-1))
  ##
  rbPal <- colorRampPalette(c('black','red'))
  cuts=length(table(links[,c(col_weight),drop=T]))
  if (cuts>1){
    #cuts=length(table(links[,c(col_weight),drop=T]))
    links$col=rbPal(cuts)[as.numeric(cut(-log10(links[,c(col_weight),drop=T]),breaks = cuts))] # coloured edges
  }
  else{
    links$col='#000000' #rbPal(1)[rep(1,nrow(links))] # all black
  }
  MyClickScript <- ifelse(isFALSE(click),NULL,
    sprintf('alert(d.name + " statistics: node size based on %s of " + d.nodesize + "; if varying, edge weights and colors represent the size of supporting literature for a given co-occurrence");',nodesize)
    )
  #
  nodes[,nodesize]=sapply(nodes[,nodesize,drop=T],function(v)(min(v,36)))
  #print(head(links))
  #print(head(nodes)) 
  #
  p=forceNetwork(Links = links,
                 Nodes = nodes,
                 Source = "V1",
                 Target = "V2",
                 NodeID = nodeid,
                 Group = group,
                 Value = edge_weight,
                 linkWidth = "0.3",#JS("function(d) { return Math.pow(d.value,1/3); }"),
                 zoom = T,
                 Nodesize = nodesize,linkDistance = JS("function(d){return 2-Math.exp(1-d.value)}"),
                 opacity = 1, legend = TRUE, opacityNoHover = 0.3,
                 bounded = F,
                 radiusCalculation = "2*d.nodesize+2", #"Math.pow(d.nodesize,1/2)+2",
                 clickAction=MyClickScript,
                 colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"),
                 linkColour=links$col)
  return(p)
}
###
if ((length(intersect(c('edge_weight','R1C1'),colnames(edges)))==2) & (length(intersect(c('category','qscore'),colnames(nodes)))==2)){
  print("Generate interactive Gene-MeSH Network based on current statistics...")
  p=widnet(
  links=edges,
  nodes=nodes,
  source=colnames(edges)[1],
  target=colnames(edges)[2],
  nodeid =colnames(nodes)[1],
  edge_weight='inv_edge_weight',
  col_weight = 'inv_edge_weight',
  group='category',
  nodesize = 'qscore')
  saveWidget(p, file=paste(wd,tag,"_interactive_Gene-MeSH_Network.html",sep=''))
}
###
NetDist=function(edges,nodes){
  print("generate graph")
  g.tm <- graph_from_data_frame(edges,directed = FALSE,vertices=nodes) # conf as edge attribute
  #CC=components(g.tm)
  # be sure it's a simple graph #
  g.tm=simplify(g.tm, edge.attr.comb="median")
  # use only 1 completely connected graph
  #CompNodes = which(CC$membership == 1)
  #g.tm = induced_subgraph(g.tm,CompNodes)
  #
  E(g.tm)$inv_edge_weight[E(g.tm)$inv_edge_weight == 0] <- 1e-15 
  #
  print(c("graph size:",gsize(g.tm)))
  ### DISTANCE MATRICES FOR KNET
  print("compute distance metrics...")
  D <- DistGraph(g.tm,edge.attr = "inv_edge_weight", dist.method = "shortest.paths", verbose=T) 
  B <- BinGraph(D,verbose=T,nsteps = 1000)
  return(list('graph'=g.tm,'B'=B, 'nodes'=V(g.tm)))
}
#
NetS=list("NEIGHnet"=NetDist(edges,nodes))
### compute Poisson occurrence for selfe loops ###
V(NetS$NEIGHnet$graph)$poiss=(1-ppois(V(NetS$NEIGHnet$graph)$occurrence,mean(V(NetS$NEIGHnet$graph)$occurrence)))/(1-ppois(0,mean(V(NetS$NEIGHnet$graph)$occurrence)))
###
print("get cliques...")
cliqueset <- max_cliques(NetS$NEIGHnet$graph,min = 3,subset = sort(V(NetS$NEIGHnet$graph)$name))
sprintf("...done, the graph contains %i cliques", length(cliqueset))
#
if (length(cliqueset)<1){
  quit(save="no",status=12)
} else{
  print("check passed: at least 1 clique")
}
#
lpaths.in.graph=lapply(cliqueset,function(g)(g$name)) #split(cliqueset$KW, cliqueset$cliqueId)
names(lpaths.in.graph)=seq(1,length(lpaths.in.graph))      
#
get.chunks=function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
#
#ncores=min(length(lpaths.in.graph),ncores)
#
#### CUSTOM LOW-MEMORY KNET FUNCTION #### 
customKnet <- function(
    g, 
    nperm=100, 
    dist.method=c("shortest.paths", "diffusion", "mfpt"), 
    vertex.attr="pheno",
    edge.attr=NULL,
    correct.factor=1,
    nsteps=1000, 
    prob=c(0, 0.05, 0.5, 0.95, 1),
    B=NULL,
    verbose=TRUE
  ){
    # calculate the Knet function for a graph, along with permutations if required

    # check that vertex weights and edge distances are present and suitable. Convert if neccessary
    dist.method <- match.arg(dist.method)
    g <- CheckAttributes(g, vertex.attr, edge.attr, verbose)
    
    # the results are saved in a list with an entry for each vertex.attr. Done even if only 1 vertex.attr supplied
    tmp <- vector("list", 2)
    names(tmp) <- c("AUK.obs","pval")
    class(tmp) <- "Knet"
    res <- rep(list(tmp), length(vertex.attr)) 
    names(res) <- vertex.attr
    nvertices <- as.integer(vcount(g))
    
    if (is.null(B)) {
        # compute B and D
        D <- DistGraph(g, edge.attr=edge.attr, dist.method=dist.method, correct.inf=TRUE, correct.factor=correct.factor, verbose=verbose) 
        B <- BinGraph(D, nsteps=nsteps, verbose=verbose) 
        rm(D)
    } else {
        # check that B is of the correct dimensions
        if (!identical(dim(B), rep(nvertices, 2))) stop("B is not of the correct dimensions")
    }
    
    # convert B to a vector and compute the maximum of B
    Bv <- as.integer(as.vector(B))
    maxB <- as.integer(max(B))
    rm(B)
    
    # run the function for each vertex attribute in vertex.attr
    for (attr in vertex.attr) {
        attr.message <- if (length(vertex.attr) == 1) NULL else paste("(", which(vertex.attr == attr), "/", length(vertex.attr), ") ", sep="")
        if (verbose) message("computing the clustering of the '", attr, "' weights ", attr.message, "using ", nperm, " permutations... ", appendLF=FALSE)
        
        # extract vertex weights 
        vertex.weights <- get.vertex.attribute(g, attr)
        vertex.weights.is.na <- is.na(vertex.weights) # don't permute the NA values
        vertex.weights[vertex.weights.is.na] <- 0
        vertex.weights <- as.double(vertex.weights)
                
        # calculate the observed netK and netAUK
        K.obs <- .Call("computenetK_fewzeros", Bv, vertex.weights, nvertices, maxB)
        res[[attr]]$AUK.obs <- sum(K.obs) / length(K.obs)
        rm(K.obs)
        # if specified, run the Knet function on permutations of the graph
        if (!is.null(nperm) & nperm > 0) {
          # run permutations - SET SEED
            K.perm <- list()
            for (i in seq_len(nperm)) {
                set.seed(i)
                vertex.weights.shuffled <- Shuffle(vertex.weights, ignore=vertex.weights.is.na)
                K.perm[[i]] <- .Call("computenetK_fewzeros", Bv, vertex.weights.shuffled, nvertices, maxB)
            }
            K.perm <- do.call("cbind", K.perm)

            # calculate the quantiles, AUK and p-values (through the z-score) for the permutations
            #res[[attr]]$K.quan <- apply(res[[attr]]$K.perm, 1, function(x) quantile(x, prob=prob))
            AUK.perm <- apply(K.perm, 2, function(x) sum(x) / length(x))
            res[[attr]]$pval <- pnorm((res[[attr]]$AUK.obs - mean(AUK.perm)) / sd(AUK.perm), lower.tail=FALSE) 
            rm(K.perm,AUK.perm)
        } else {
            # if no permutations are run, permutation-related statistics are returns equal to NA
            res[[attr]][c("pval")] <- NA
        }
         
        if (verbose) message(" done")
    }
    rm(Bv)
    gc(verbose=F)
    # cleanup and output
    #if (length(vertex.attr) == 1) res <- res[[1]] # if only one vertex attribute is input, don't return a list of lists of results
    return(res)
}
####
environment(customKnet) <- asNamespace('SANTA')
assignInNamespace("Knet", customKnet, ns = "SANTA")
####
Kstat=function(chunk){
  pws=names(chunk)
  badp=c()
  exp.nodes=V(g.tm)$name
  for (pw in pws) {
    if (sum(as.numeric(unname(chunk[pw][[1]]) %in% exp.nodes)) >= 3){
      g.tm <- set_vertex_attr(
        g.tm, name=pw,
        value=as.numeric(exp.nodes %in% unname(chunk[pw][[1]]))) 
    }
    else{
      badp=c(badp,pw)
    }
  }
  #
  pws=setdiff(pws,badp)
  
  
  ###
  #set.seed(2202)
  result_knet <- SANTA::Knet(g.tm, nperm=500, vertex.attr=pws, edge.attr=NULL, verbose=F,B = B)
  #result_knet <- Compactness(g.tm, nperm=300, vertex.attr=pws, edge.attr="weight", verbose=T)
  # if (length(pws) == 1){
  #   result_knet=list(result_knet)
  #   names(result_knet)=sprintf("%s",pws[1])
  # }
  rm(g.tm,B,pos = environment())
  gc(verbose = F)
  return(result_knet)
}
#
for (data in names(NetS)){
  print("estimate available memory...")
  gc(verbose=F)
  free=as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern=TRUE))
  #free=as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE)) # available memory at that time
  use=0.25*free
  print(sprintf("set to use around 25%% of the available RAM, i.e. %g GB", round(use/1E6)))  
  net=NetS[data][[1]]
  g.tm=net$graph
  B=net$B
  # object.size() returns bytes
  # needs a 2x moltiplication factor because R parallel sucks!
  # each core will need
  # 1) A bin matrix B
  # 2) The graph object (g.tm)
  # 3) a vector of size max(B) that contains the observed Knet values for each bin (to compute AUK)
  # 4) a matrix nperm X max(B) that contains sampled Knet values 
  # 5) an nper vector that containst sampled AUK values  
  # 6) at most each core needs to evaluate #cliques objects 
  required=as.numeric(object.size(B)+object.size(g.tm)+
  object.size(numeric(500))+object.size(numeric(max(NetS$NEIGHnet$B)))+object.size(matrix(1000,500))+
  object.size(rep(vector('list',2),length(lpaths.in.graph))))/1E3 # returns KB
  print(sprintf("given B and 500 Knet permutations,each core will need %f GB",required/1E6))
  # set how many cores will then be used
  mcores=max(min(floor(use/required),ncores),1)
  print(sprintf("Then %g cores will be used",mcores))
  print("initialize parallelization")
  if (mcores>1){
    chunks=get.chunks(lpaths.in.graph,min(length(lpaths.in.graph),mcores))
    cl <- makeCluster(length(chunks),type = 'PSOCK')
    clusterEvalQ(cl, library(SANTA))
    clusterExport(cl, c("g.tm","B","Kstat","customKnet"))
    clusterEvalQ(cl,environment(customKnet) <- asNamespace('SANTA'))
    clusterEvalQ(cl,assignInNamespace("Knet", customKnet, ns = "SANTA"))
    #
    print(c("Knet applied to",data))
    result_knet=parLapplyLB(cl, chunks, Kstat)
    stopCluster(cl)
  }else{
    chunks=list('1'=lpaths.in.graph)
    result_knet=lapply(chunks, Kstat)   
  }
  print(c("Knet completed",data))
  print("Prepare results table...")  
  #print(str(result_knet))
  result_knet=do.call("c", result_knet)
  #print(result_knet)
  names(result_knet)=str_replace_all(names(result_knet),"[0-9]+\\.","")  
  #pvals=sapply(result_knet, function(res){return(res$pval)})
  pvals=sapply(result_knet,'[[','pval')# function(res){return(res$pval)})
  #pvals=sort(pvals, decreasing = FALSE)
  AUK.obs=sapply(result_knet,'[[','AUK.obs')
  enrichment<- data.frame(names(pvals), unname(pvals), p.adjust(unname(pvals), method = 'BH'),AUK.obs, row.names = NULL,stringsAsFactors = F)
  colnames(enrichment)=c('cliqueID','pKnet','pAdjust','AUK.obs')
  rm(result_knet)
  gc(verbose=F)
  print("...done")
}  
#
### select significant cliques ###
hit.cliques=lpaths.in.graph[enrichment[enrichment$pAdjust<=0.01,]$cliqueID]
# be sure of sorting hit.cliques by name #
hit.cliques=hit.cliques[as.character(sort(as.numeric(names(hit.cliques))))]
hit.cliques.AUK=setNames(enrichment[enrichment$pAdjust<=0.01,]$AUK.obs,names(hit.cliques))
### record Knet 
#
if (length(hit.cliques)<2){
  quit(save="no",status=14)
} else{
  print("check passed: at least 2 significant clique")
}
#
hot.KW=sort(table(do.call("c",hit.cliques)),decreasing = T)
#
### check shift in average clique size after Knet statistics ###
data=reshape2::melt(list(`post-Knet` = sapply(hit.cliques, length), `pre-Knet` = sapply(lpaths.in.graph, length)))
colnames(data)=c('value','checkpoint')
ggplot(data,aes(x=value, fill=checkpoint)) + geom_density(alpha=0.2) + ggtitle("Shift effect of Knet statistics on the size of significant cliques") + labs(y='density',x="clique size")
ggsave(paste(wd,tag,"_Knet_shift_effect.png",sep=''),width = 12, height = 12,dpi = 320)
###
### PRUNE NETWORK BY HOT KEYWORDS, THEN COMMUNITY DETECTION ### 
g=NetS$NEIGHnet$graph
### USE SIMPSON'S SIMILARITY AND STOUFFER'S COMBINED P-VALUE METHOD ###
print("compute co-occurrence based probability score and use it as edge weight...")
g <- set.edge.attribute(g, "weight", index=E(g), as.numeric(format(1-E(g)$inv_edge_weight,digits=15))) #(E(g)$R1C1-mean(E(g)$R1C1))/sd(E(g)$R1C1))
#minZ=(0-mean(E(g)$R1C1))/sd(E(g)$R1C1)
#maxZ=max(E(g)$weight)#(1-mean(E(g)$simpson))/sd(E(g)$simpson)
#dicesim=similarity(g,method='dice')
#colnames(dicesim)=V(g)$name;rownames(dicesim)=V(g)$name
# Simpson = Dice * sum(degs)/min(deg1,deg2) #
#cl <- makeCluster(ncores,type = 'PSOCK')
#clusterExport(cl, c("g"))
#clusterEvalQ(cl,library(igraph))
#E(g)$simpson=dicesim[ends(g,E(g))]*parApply(cl,ends(g,E(g)),1,function(v)(sum(unname(degree(g,v)))))/(2*parApply(cl,ends(g,E(g)),1,function(v)(min(unname(degree(g,v))))))
#E(g)$dice=dicesim[ends(g,E(g))]
#E(g)$invlog=dicesim[ends(g,E(g))]
# infer correlation between nodes from distances #
print("...and use it to infer distance-based, node-node correlations")
#d.tm=distances(g,V(g)$name,V(g)$name,weight=E(g)$simpson)
#d.tm[!(is.finite(d.tm))]=max(d.tm[is.finite(d.tm)])+1
#clusterExport(cl, c("d.tm"))
#cor.tm=parLapply(cl,combn(V(g)$name,2,simplify=F),function(l){
#  d1=d.tm[l[[1]],,drop=T];d2=d.tm[l[[2]],,drop=T]
  # restrict to common neighbors #
#  oknames=intersect(names(d1[is.finite(d1)]),names(d2[is.finite(d2)]))
#  return(cor(unname(d1[match(oknames,names(d1))]),unname(d2[match(oknames,names(d2))])))
#})
#
# symmdist=function(objs,distcom,fillnav,filld=0){
#   mdist=matrix(nrow=length(objs),ncol = length(objs))
#   colnames(mdist)=objs;rownames(mdist)=objs  
#   mdist[lower.tri(mdist)]=distcom
#   mdist=t(mdist)
#   mdist[lower.tri(mdist)]=distcom
#   diag(mdist)=filld
#   mdist[is.na(mdist)]=fillnav
#   return(mdist)
# }
# #
# cor.tm=symmdist(V(g)$name,unlist(cor.tm),fillnav=0,filld=1)
# ###
# cor.graph=graph_from_adjacency_matrix(matrix(1,nrow=nrow(cor.tm),ncol=ncol(cor.tm),dimnames=list(rownames(cor.tm),colnames(cor.tm))))
# E(cor.graph)$correlation=apply(ends(cor.graph,E(cor.graph)),1,function(v)(cor.tm[v[1],v[2]]))
###
#plot(g,layout=layout_nicely,vertex.size=4,vertex.label=NA)
g=induced_subgraph(g,names(hot.KW))
#
### "x" consists of a collapsed distance measure between members of two cliques. Used to build a clique-clique graph ###
#plot(g,layout=layout_nicely,vertex.size=4,vertex.label=NA)
# EDGE WEIGHT HAS TO BE RE-SET!
# following the method from "Assigning protein function from domain‑function associations using DomFun"
# we convert simpsons score to Z-scores then we will use their associated pvalues
### compute maxZ and minZ for the current mean and sd ###
###
# using geometric mean 
# using indicator function 
# print("calculate intersections of cliques...")
# clusterExport(cl, c("g"))
# ### SMALLER MEANS CLOSER ###
# g.d=distances(g)#,weight=as.numeric(format(1-E(g)$inv_edge_weight,digits=15)))
# g.a=as_adj(g)
# g.da=as.matrix(g.d*g.a)
# diag(g.da)=1#-setNames(V(g)$poiss,V(g)$name)[colnames(g.da)]
# #
# scores=function(R,C,d.neigh=g.da){
#   dij=d.neigh[R,C]
#   # if the network is not complete, infinite distances would cause NaN values
#   if (dij==0 | is.nan(dij)){
#     return(1e-15)
#   }else{
#     return(dij)
#     }}
# #
# score_v=Vectorize(scores)
# #

# #
# clusterExport(cl,c('score_v','g.da'))#,'maxZ','minZ'))
# #
# com <- vector(mode = "list", length = choose(length(hit.cliques),2))
# com <-combn(hit.cliques,2,simplify=F)
# x=vector(mode = "list", length = length(com))
# print("...by combining Fisher-Pearson combined pvalue statistics - see Heard 2017...")
# set.seed(2202)
# ## code for testing 
# l=com[names(sapply(com,'[',1))=='2' & names(sapply(com,'[',2))=='3']
# ##
# x=parLapply(cl,com,function(l){
#   #
#   print(l)
#   #m=outer(setNames(nm=l[[1]]),setNames(nm=l[[2]]),score_v)
#   m=g.da[match(l[[1]],rownames(g.da)),match(l[[2]],rownames(g.da))]
#   # # include correlation between test statistics #
#   #vec <- unique(c(l[[1]],l[[2]]))
#   #cor.graph=graph_from_adjacency_matrix(matrix(1,nrow=length(vec),ncol=length(vec),dimnames=list(vec,vec)))
#   #E(cor.graph)$correlation=apply(ends(cor.graph,E(cor.graph)),1,function(v)(cor.tm[v[1],v[2]]))
#   # test1=outer(l[[1]],l[[2]],function(x,y) paste(x,y,sep='_'))
#   # test2=outer(as.vector(test1),as.vector(test1),Vectorize(function(x,y){
#   #   l1=stringr::str_split(x,'_')[[1]]
#   #   l2=stringr::str_split(y,'_')[[1]]
#   #   inter=intersect(l1,l2)
#   #   #return(mean(E(cor.graph)[setdiff(l1,inter) %--% setdiff(l2,inter)]$correlation))
#   #   return(length(inter)/length(unique(c(l1,l2))))
#   #   }))
#   # diag(test2)=1
#   # test2[!(is.finite(test2))]=0
#   # # avoid 1/0s #
#   #ps=pnorm(as.vector(m[!(is.na(m))]))
#   ps=as.vector(m)
#   ps=ifelse(ps==0 | is.nan(ps),1e-15,ps)
#   ps=ifelse(ps==1,as.numeric(format(1-1e-15,digits=15)),ps)  
#   # combined Fisher/Pearson with w=0.5=(1-a)/(b-a), see Heard 2017 #
#   #SfSp=2*(sum(log(ps))-sum(log(1-ps)))
#   #SfSp=function(ps)(list(statistic=(0.5*sum(log(ps))-0.5*sum(log(1-ps))),p=pnorm(0.5*sum(log(ps))-0.5*sum(log(1-ps)))))
#   #nulldist=qchisq(nullps,df=2*length(ps))-qchisq(1-nullps,df=2*length(ps))
#   stouf=poolr::fisher(1-ps)#,adjust='generalized',R=poolr::mvnconv(test2))#,target='z',side=1,cov2cor=T)) 
#   #stoufgaoz=poolr::stouffer(pnorm(as.vector(m[!(is.na(m))])),adjust='gao',R=poolr::mvnconv(test2,target='z',side=1,cov2cor=T))   
#   #df=data.frame(names(l)[1],names(l)[2],stouf$p,stouf$statistic,stringsAsFactors=F)
#   # scores=mapply(function(g1,g2)(paste(g1,g2)),cross[,1,drop=T],cross[,2,drop=T])
#   # #
#   #inter <- intersect(l[[1]],l[[2]])
#   #vec <- unique(l[[1]],l[[2]]) 
#   #Ivec <- 0+(vec %in% inter)#(l[[1]] %in% inter)+(l[[2]] %in% inter) #-length(inter) # indicator function
#   #df=data.frame(names(l)[1],names(l)[2],mean(c(igraph::E(g)[l[[1]] %--% l[[2]]]$inv_edge_weight,Ivec)),stouf$p,stouf$statistic,stringsAsFactors=F)
#   df=data.frame('c1'=names(l)[1],'c2'=names(l)[2],'stouf.p'=stouf$p,'weight'=stouf$statistic,stringsAsFactors=F)
#   #colnames(df)=c('c1','c2','stouf.p','weight')
#   return(df)
# })
# #
# stopCluster(cl)
# rm(com)
# gc(verbose=F)
# print("...done")
# ### immediately remove edges with 0 weight ###
# #x=x[sapply(x,function(df)(df$weight!=0))]
# #
# n=length(x)
# dtab <- data.table(c1=rep(NA,n), c2=rep(NA,n),weight=rep(0,n))
# dtab=rbindlist(x)
# dtab$p.adjust=p.adjust(dtab$stouf.p,'BH')
# dtab.keep=dtab[dtab$p.adjust<=0.01]
# #colnames(x)=c('c1','c2','weight')
# #wecdf=Vectorize(ecdf(x$weight))
# #x$weight=-log10(1-x$weight)
# #
# ### CLIQUE-CLIQUE GRAPH ###
# print("generate cliques network...")
# GC=graph_from_data_frame(dtab.keep,directed = F)
# print("clean environment...")
# rm(x,dtab,pos=environment())
# ###
# print(sprintf("%i cliques were pruned, as no significant clique-clique edge was found", length(setdiff(names(hit.cliques),V(GC)$name)
# )))
# #
# hit.cliques=hit.cliques[names(hit.cliques) %in% V(GC)$name]
# hit.cliques.AUK=hit.cliques.AUK[names(hit.cliques.AUK) %in% V(GC)$name]
# ### not all cliques may have a significant connection to other cliques ### 
# ##
# gc(verbose=F)
# ### PRUNE EDGES WITH 0 WEIGHT ###
# #print("generate cliques network...")
# #GC=delete.edges(GC, which(E(GC)$weight < quantile(E(GC)$weight,0.5)))
# #CC=cluster_louvain(GC) 
# CC=cluster_leiden(GC,objective_function = 'modularity',beta=0,n_iterations=10000)# ,weights = NULL,resolution_parameter = 1.5,beta = 0.001,n_iterations = 1000)
# #
# #CC=cluster_leiden(g,objective_function = 'modularity',n_iterations = 1000)
# #
# if (length(communities(CC)) < 2){
#   quit(save="no", status=13)
# } else{
#   print("check passed: at least 2 communities")
# }
# #
# print(sprintf("the graph's modularity is %f",modularity(GC,weight=E(GC)$weight,membership = CC$membership)))
# png(paste(tag,"_clique-clique_communities_net.png",sep=''),width = 12, height = 12,units = 'in',bg='grey95',res=200)
# plot(GC, vertex.color=rainbow(max(CC$membership), alpha=0.6)[CC$membership],layout=layout.graphopt,edge.width=log10(E(GC)$weight),vertex.label.cex=0.65)+title("communities in cliques-only graph")
# dev.off()
# ### Propagate memberships to individual entities ### 
# clust_named=communities(CC)
# ###
# print("check Cstack Info...")
# print(Cstack_info())
# ### sugiyama ###
# layex <-  layout_with_sugiyama(GC, layers=apply(sapply(clust_named,
#                                                        function(x) V(GC)$name %in% as.character(x)),
#                                                 1, which))
# #
# origvert <- c(rep(TRUE, vcount(GC)), rep(FALSE, nrow(layex$layout.dummy)))
# realedge <- as_edgelist(layex$extd_graph)[,2] <= vcount(GC)
# png(paste(tag,"_clique-clique_communities_Sugiyama.png",sep=''),width = 16, height = 12,units = 'in',bg='grey95',res=200)
# plot(layex$extd_graph,
#      edge.arrow.size=.1,
#      edge.width=log10(E(GC)$weight),vertex.label.cex=0.65,
#      vertex.color=rainbow(max(CC$membership), alpha=0.7)[CC$membership],
#      vertex.size=5,
#      vertex.shape='circle',
#      vertex.label=NA,
#      sub='Use it to visualize connectivity between communities'
# )+title("cliques-cliques Sugiyama graph")
# dev.off()
# rm(layex,pos=environment())
# gc(verbose=F)
#### CLUSTER PROPAGATION #### 
##### SET SEED ON LEIDEN CLUSTERING ######
custom_leiden=function (graph, objective_function = c("CPM", "modularity"), 
  weights = NULL, resolution_parameter = 1, beta = 0.01, initial_membership = NULL, 
  n_iterations = 2, vertex_weights = NULL) 
{
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  objective_function <- igraph.match.arg(objective_function)
  objective_function <- switch(objective_function, cpm = 0, 
    modularity = 1)
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && !any(is.na(weights))) {
    weights <- as.numeric(weights)
  }
  else {
    weights <- NULL
  }
  if (!is.null(initial_membership) && !any(is.na(initial_membership))) {
    initial_membership <- as.numeric(initial_membership)
  }
  else {
    initial_membership <- NULL
  }
  if (!is.null(vertex_weights) && !any(is.na(vertex_weights))) {
    vertex_weights <- as.numeric(vertex_weights)
    if (objective_function == 1) {
      warning("Providing node weights contradicts using modularity")
    }
  }
  else {
    if (objective_function == 1) {
      vertex_weights <- strength(graph, weights = weights)
      resolution_parameter <- resolution_parameter/sum(vertex_weights)
    }
  }
  on.exit(.Call(C_R_igraph_finalizer))
  membership <- initial_membership
  if (n_iterations > 0) {
    for (i in 1:n_iterations) {
      set.seed(i)
      res <- .Call(C_R_igraph_community_leiden, graph, 
        weights, vertex_weights, as.numeric(resolution_parameter), 
        as.numeric(beta), !is.null(membership), membership)
      membership <- res$membership
    }
  }
  else {
    prev_quality <- -Inf
    quality <- 0
    while (prev_quality < quality) {
      prev_quality <- quality
      res <- .Call(C_R_igraph_community_leiden, graph, 
        weights, vertex_weights, as.numeric(resolution_parameter), 
        as.numeric(beta), !is.null(membership), membership)
      membership <- res$membership
      quality <- res$quality
    }
  }
  res$algorithm <- "leiden"
  res$vcount <- vcount(graph)
  res$membership <- res$membership + 1
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    res$names <- vertex_attr(graph, "name")
  }
  class(res) <- "communities"
  res
}
environment(custom_leiden) <- asNamespace("igraph")
#assignInNamespace("cluster_leiden", custom_leiden, ns = "igraph")    
#####
i=0
set.seed(i);tCC=cluster_leiden(g,objective_function = 'modularity',n_iterations=1)# ,weights = NULL,resolution_parameter = 1.5,beta = 0.001,n_iterations = 1000)
#tCC=cluster_leiden(g,resolution_parameter=1/(2*sum(strength(g))),vertex_weights=strength(g),beta=0,n_iterations=1)
newmod=modularity(g,weight=E(g)$weight,membership = tCC$membership)
lastmod=0
while ((newmod-lastmod)>0){
  print(sprintf("last best modularity: %f",newmod))
  if (exists('newtCC')){
    tCC=newtCC
  }
  i=i+1
  set.seed(i);newtCC=cluster_leiden(g,objective_function = 'modularity',n_iterations=1,initial_membership=tCC$membership)
  lastmod=newmod
  newmod=modularity(g,weight=E(g)$weight,membership = newtCC$membership)  
}
print(sprintf("found local maximum in %i iterations",i))
# for (i in seq(1E2)){
#   set.seed(2202)
#   tCC=cluster_leiden(g,objective_function = 'modularity',beta=0,n_iterations=1,initial_membership=tCC$membership)# ,weights = NULL,resolution_parameter = 1.5,beta = 0.001,n_iterations = 1000)  
#   #tCC=cluster_leiden(g,resolution_parameter=1/(2*sum(strength(g))),vertex_weights=strength(g),beta=0,n_iterations=1,initial_membership=tCC$membership)
# }
#tCC=cluster_leiden(g,resolution_parameter=1/(2*sum(strength(g))),vertex_weights=strength(g),beta=0,n_iterations=10000)# ,weights = NULL,resolution_parameter = 1.5,beta = 0.001,n_iterations = 1000)
print(sprintf("the graph's modularity is %f",modularity(g,weight=E(g)$weight,membership = tCC$membership)))
#
clust_prop=communities(tCC)
#lapply(clust_named,function(C)(unique(do.call('c',hit.cliques[C]))))
for (i in names(clust_prop)){
  clust_prop[[i]]=data.frame('entity'=clust_prop[[i]],'membership'=i,stringsAsFactors = F)
}
clust_prop=bind_rows(clust_prop)
clust_prop$membership=as.numeric(clust_prop$membership)#+1
# MULTIMAPPING WILL HAVE THEIR OWN MEMBERSHIP
print("plot community size distribution...")
multi=as.character(unique(clust_prop[duplicated(clust_prop$entity),]$entity))
sizes=table(clust_prop$membership)
tabs=table(sizes)
png(filename=paste(tag,"_cliques_community_size_distribution.png",sep=''),width = 12, height = 8,units = 'in',bg='white',res=200)#,bg='transparent')
#barplot(sort(sizes,decreasing = T),space=.2,col = rainbow(max(as.numeric(names(sizes))), alpha=0.7),xlab='Community ID',ylab='Count',main='Cliques Community Size Distribution')
barplot(sort(sizes,decreasing = T),space=.2,col = 'gray80',ylim=c(0,max(sizes)),cex.lab=2.5,cex.main=2.5,xpd=F,cex.axis=2.5,,xlab='Community ID',ylab='Count',main=sprintf('Cliques Community Size Distribution\nModularity=%.2f',modularity(g,weight=E(g)$weight,membership = tCC$membership)))
dev.off()
#print(sprintf("there are %g ties among community sizes", sum(tabs[tabs>1])))
#
#memberships=tapply(clust_prop$membership,factor(clust_prop$entity),function(v)(v))
#memberships=ifelse(sapply(memberships,function(v)(length(v)>1)),'multimapping',sapply(memberships,'[[',1))
#clust_prop=data.frame('entity'=names(memberships),'membership'=unname(memberships))
#print(tabs)
#print(multi)
## iterative assignment ###
# for (i in multi){
#   #print(i)
#   #set.seed(2202)
#   cbymin=sort(names(sizes[sizes==min(sizes[names(sizes) %in% clust_prop[clust_prop$entity==i,]$membership])]))
#   if (length(cbymin)>1){
#     # then decide by AUK.obs and if not by name #
#     top.clique=names(which.max(hit.cliques.AUK[names(hit.cliques)[sapply(hit.cliques,function(hc)(any(i %in% hc)))]]))
#     comm2belong=names(clust_named[sapply(clust_named,function(C)(top.clique %in% C))])
#     clust_prop[clust_prop$entity == i,]$membership=comm2belong  
#   }else{
#     clust_prop[clust_prop$entity == i,]$membership=cbymin  
#   }
#   #
#   sizes=table(clust_prop$membership)
#   #print(sizes)
# }
#clust_prop[clust_prop$entity %in% multi,]$membership=sapply(clust_prop[clust_prop$entity %in% multi,]$entity,
#                                                            function(e){
#
clust_prop=clust_prop[!(duplicated(clust_prop)),]
cnodes=clust_prop
cnodes$nodesize=1
cnodes=cnodes[!(duplicated(cnodes)),]
clust_prop_df=clust_prop                                                              
### RE-SUBSET co-occ net
g=NetS$NEIGHnet$graph
knodes_named=setNames(nodes$qscore,nodes$KW)
g=set_vertex_attr(g,"qscore", value = knodes_named[match(V(g)$name,names(knodes_named))])
g=induced_subgraph(g,cnodes$entity)
clust_prop=setNames(as.character(clust_prop$membership),clust_prop$entity)
### plot clique memberships propagated on entities
png(paste(tag,"_gene-MeSH_CliqueNetwork_by_clique_communities.png",sep=''),width = 12, height = 12,units = 'in',bg='grey95',res=200)
plot(g, vertex.color=ifelse(unname(clust_prop[V(g)$name])=='multimapping','grey93',rainbow(max(as.numeric(as.factor(clust_prop))), alpha=0.6)[as.numeric(unname(clust_prop[V(g)$name]))]),layout=layout.graphopt,niter=500,weigths=-log(E(g)$inv_edge_weight),niter=500,edge.width=5*E(g)$weight,vertex.label.cex=0.5,sub='ambiguous nodes are assigned to the smallest communities and highest-ranking cliques')+title("Gene-MeSH, cliques-only network coloured by clique-based communities")
dev.off()
###

#### SAVE IT ALSO AS INTERACTIVE FORMAT #### 
print("generate interactive network of cliques...")
### convert igraph ###
V(g)$group=unname(clust_prop[V(g)$name])
g_toD3=igraph_to_networkD3(g, group = V(g)$group)
### transfer edge attributes ###
cl <- makeCluster(ncores,type = 'PSOCK')
clusterExport(cl,c('g','g_toD3'))#,'maxZ','minZ'))
clusterEvalQ(cl,library(igraph))
###
g_toD3$links$edge_weight=parApply(cl,g_toD3$links[,c(1,2)],1,function(v){
  v=as.character(g_toD3$nodes[v+1,]$name)
  return(E(g)[v[1] %--% v[2]]$inv_edge_weight)
  })
g_toD3$links$inv_edge_weight=parApply(cl,g_toD3$links[,c(1,2)],1,function(v){
  v=as.character(g_toD3$nodes[v+1,]$name)
  return(E(g)[v[1] %--% v[2]]$inv_edge_weight)
  })
#
g_toD3$links$R1C1=parApply(cl,g_toD3$links[,c(1,2)],1,function(v){
  v=as.character(g_toD3$nodes[v+1,]$name)
  return(E(g)[v[1] %--% v[2]]$R1C1)
  })
# keep nodesize to 1
g_toD3$nodes$nodesize=knodes_named[match(g_toD3$nodes$name,names(knodes_named))]
#
rbPal <- colorRampPalette(c('black','red'))
cuts=length(table(g_toD3$links$inv_edge_weight))
if (cuts>1){
  #cuts=length(table(links[,c(col_weight),drop=T]))
  g_toD3$links$col=rbPal(cuts)[as.numeric(cut(-log10(g_toD3$links$inv_edge_weight),breaks = cuts))] # coloured edges
}else{
  g_toD3$links$col='#000000'}
#
MyClickScript <- sprintf('alert(d.name + " statistics: node size based on %s of " + d.nodesize + "; if varying, edge weights and colors represent the size of supporting literature for a given co-occurrence");','nodesize')
#
### COLOR PALETTE ACCORDING TO NUMBER OF DIFFERENT COMMUNITIES ###
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
print('#### color vector ####')
print(col_vector[seq(length(unique(V(g)$group)))])
colpal=JS(paste0('d3.scaleOrdinal(',paste0('[',paste(shQuote(col_vector[seq(length(unique(V(g)$group)))], "cmd"),collapse=','),']'),');'))
#
p=forceNetwork(
  Links=g_toD3$links,
  Nodes=g_toD3$nodes,
  Source='source',
  Target='target',
  NodeID ='name',
  Value='edge_weight',
  linkWidth = "0.3",#JS("function(d) { return Math.pow(d.value,1/3); }"),
  linkColour = g_toD3$links$col,
  Group='group',
  Nodesize = 'nodesize',
  zoom=T,
  opacity = 1, legend = TRUE, opacityNoHover = 0.3,
  bounded = F,
  radiusCalculation = "2*d.nodesize+2", #"Math.pow(d.nodesize,1/2)+2",
  colourScale = colpal,clickAction=MyClickScript)
#
saveWidget(p, file=paste(wd,tag,"_interactive_Cliques_Network.html",sep=''))
### build a function that constructs ranked queries ###
knodes_named=setNames(nodes$qscore,nodes$KW)
cknodes_named=knodes_named[names(knodes_named) %in% names(hot.KW)]
#clust_named=communities(CC)
clust_named=split(names(clust_prop),clust_prop)
#
### INSPECT cluster sizes ### 
print("check clusters with viable sizes...")
okCC=list()
mockcomb=combs
while (mockcomb >=3){
  okCC=clust_named[sapply(clust_named,function(v)(length(v)>=mockcomb))]
  if (length(okCC)<length(clust_named)){
    mockcomb=mockcomb-1
  }
  else{
    break
  }
}
#
if (length(okCC)==0){
  stop(sprintf("FAIL: NO COMMUNITY HAS A MINIMAL SIZE OF 3 ENTITIES")) 
} else if (mockcomb != combs){
  print(sprintf("WARNING: SOME COMMUNITIIES' SIZE IS SMALLER THAN %i:",combs))
  print(clust_named[sapply(clust_named,function(v)(length(v)<combs))])
} else{
  print(sprintf("ALL COMMUNITIIES HAVE SIZE EQUAL OR BIGGER THAN %i",combs))
}
#
#clust_named=okCC
### SAVE nodes-to-cliquecommunity table ###
write.table(clust_prop_df,paste(tag,"_nodes_in_clique-clique_communities_net.tsv",sep=''), sep="\t", row.names = F)
colnames(clust_prop_df)=c('Entity','CommunityMembership')
nodes=merge(nodes,clust_prop_df,by.x = 'KW',by.y = 'Entity',all.x=T)
nodes$CommunityMembership=ifelse(is.na(nodes$CommunityMembership),-1,nodes$CommunityMembership)
write.table(nodes,node_files[1], sep="\t", row.names = F)
#
rm(clust_prop,pos=environment())
gc(verbose=F)
#### THE CONQUER ALGORITHM ####
#
# 1) FIND ENTITIES THAT INTERCONNECT COMMUNITIES 
# 2) FIND THEIR NEIGHBOURS 
# 3) BASED ON KNODE AND EDGE WEIGHTS, RANK QUERIES OF SIZE COMB. THESE QUERIES WILL CONTAIN TWO INTERCONNECTING ENTITIES + COMB-2 NEIGHBOURS
#
print("apply the CONQUER algorithm to bridge communities by querying...")
#
### LOG TRANSFORM EDGE AND VERTEX WEIGHTS ### 
E(g)$weight=ifelse(E(g)$weight==0,-log10(1e-15),-log10(E(g)$weight))
knodes_named=ifelse(knodes_named==0,-log10(1e-15),-log10(knodes_named))
# 
print("initialize parallelization...")
#cl <- makeCluster(ncores,type = 'PSOCK')
#clusterEvalQ(cl, library(igraph))
clusterExport(cl, c("knodes_named","max_attempts","combs","g"))
#
print("CONQUER!")
#
link=combn(okCC,2,simplify = F,function(l){
  # l is a pair of communities
  # le are the vertices that share an edge with each others' community members
  #print(l)
  le=attr(E(g)[l[[1]] %--% l[[2]]],'vnames')
  if (length(le) > 0){
    #
    lv=str_split(le,'\\|')
    #gg=induced_subgraph(g,c(l[[1]],l[[2]]))
    clusterExport(cl, c("l"),envir=environment())
    #
    print(sprintf("evaluating %i keywords-seeds linking communities %s...",length(lv),paste(names(l),collapse = ' and ')))
    # look at the neighbourhood of each seed
    ln=names(l)#paste(names(l),collapse = ' and ')
    #
    #CC_queues=parLapply(cl,X = lv,function(v){
    CC_queues=clusterApplyLB(cl,x = lv,fun=function(v){
      #print(v)
      gg=induced_subgraph(g,c(l[[1]],l[[2]]))
      nlv=V(gg)$name[neighbors(gg,v)]
      nlv=setdiff(nlv,v)
      # such neighbourhood must only contain members of the two communities
      nlv=nlv[nlv %in% c(l[[1]],l[[2]])]
      if (length(nlv)>combs-2){
        #print(nlv)
        queues=combn(nlv,combs-2,simplify=F,function(s)(c(v,unlist(s))))
        #scorele=E(g)[v[1] %--% v[2]]$weight*knodes_named[names(knodes_named) %in% v]
        #scorenle=lapply(nlv,function(n)(setNames(sum(E(g)[c(v,nlv) %--% n]$weight)*knodes_named[names(knodes_named)==n],n)))
        #rankings=combn(scorenle,combs-2,simplify=F,function(s)(c(scorele,unlist(s,use.names = T))))
        rankings=lapply(queues, function(q)(unlist(lapply(q,function(n)(setNames(sum(c(E(gg)[q %--% n]$weight,knodes_named[names(knodes_named)==n])),n))))))
        print(rankings)
        #return(attr(V(g)[v],"Knode"))
        return(rankings)  
      }
      else{
        return(list())
      }      
    })
    #
    CC_queues=sort(do.call('c',lapply(unlist(CC_queues,recursive = F),function(query)(setNames(sum(query),paste(sort(names(query)),sep='_',collapse = '_'))))),decreasing = T)
    if (!(is.null(CC_queues))){      
      CC2att=CC_queues[!(duplicated(CC_queues))][1:max_attempts]
      CC2df=data.frame(rep(ln[1],length(CC2att)),rep(ln[2],length(CC2att)),CC2att,names(CC2att),row.names = NULL,stringsAsFactors = F)
      colnames(CC2df)=c('c1','c2','weight','query')
      return(CC2df)  
    }
    else{
      CC2df=as.data.frame(matrix(ncol=4))
      colnames(CC2df)=c('c1','c2','weight','query')
      return(CC2df)  
     }
  }
  else{
    CC2df=as.data.frame(matrix(ncol=4))
    colnames(CC2df)=c('c1','c2','weight','query')
    return(CC2df)
  }
})
#
print("CONQUERED!")
stopCluster(cl)
### a network of possible queries ###
conqnet=bind_rows(link)
conqnet=conqnet[complete.cases(conqnet),]
write.table(conqnet,paste(wd,tag,"_ordered_queries.tsv",sep=''),row.names = F,col.names = T,quote = F,sep = '\t')
# #
# querank=function(commy){
#   commy.combs=combn(commy,min(length(commy),combs),simplify = F,function(v)(setNames(sum(knodes_named[names(knodes_named) %in% v]),paste(v,sep='_',collapse = '_'))))
#   return(sort(do.call('c',commy.combs),decreasing = T)[1:max_attempts])
# }
# #
# cl <- makeCluster(min(ncores,length(clust_named)),type = 'PSOCK')
# clusterExport(cl, c("knodes_named","querank","max_attempts","combs"))
# queries=parLapply(cl,clust_named,querank)
# stopCluster(cl)
# #
# queries=do.call(rbind.data.frame, lapply(queries,names))
# colnames(queries)=seq(1,ncol(queries))
# #
# write.table(queries,paste(wd,tag,"_ordered_queries.tsv",sep=''), sep="\t", row.names = F)