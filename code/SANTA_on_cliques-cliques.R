#!/home/musellla/miniconda3/envs/R363/lib/R/bin/Rscript
# SANTA RANKINGS BASED ON JACCARD DISTANCE-DERIVED NODE AND EDGE ATTRIBUTES
###
options(Ncpus = 8)
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
    nncores=floor(numargs[1])
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
list.of.packages <- c("igraph","BiocManager","tidyr","parallel","reshape2","stringr","ggplot2","snow","dplyr","networkD3","htmlwidgets","data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages(lib.loc=libloc)[,"Package"])]
if(length(new.packages)){
  options(Ncpus = 8)
  install.packages(new.packages,repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc) # FAU mirror
}
#if(!(c("SANTA") %in% installed.packages(lib.loc=libloc,force=T)[,"Package"])) BiocManager::install("SANTA",lib.loc=libloc,destdir = libloc, force=T) # FAU mirror
if(!(c("SANTA") %in% installed.packages(lib.loc=libloc,force=T)[,"Package"])){
  options(Ncpus = 8)
  install.packages("/home/musellla/txtmining/code/SANTA_2.24.0.tar.gz", repos=NULL) # FAU mirror  
} 
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
#
library(data.table)
###
#### LOAD DATA ####
if ("incomplete" %in% args){
  #print(paste(wd,"data/gene-subgraphs/",sep = '',collapse = ''))
  setwd(paste(wd,"data/",sep = '',collapse = ''))
  edge_files = unique(grep("allxall_sig_combinations_bh_to_cytoscape.tsv",list.files(),value = T)) # list.files() %>% str_subset(pattern = paste(tag,"_edges\\w+",sep = '',collapse = ''))
  node_files = unique(grep("_allxall_nodes_table_bh_to_cytoscape.tsv",list.files(),value = T)) #list.files() %>% str_subset(pattern = "nodes\\w+")
} else {
  setwd(wd)
  edge_files = unique(grep("Complete_edges_table",list.files(),value = T)) # list.files() %>% str_subset(pattern = paste(tag,"_edges\\w+",sep = '',collapse = ''))
  node_files = unique(grep("Complete_nodes_table",list.files(),value = T)) #list.files() %>% str_subset(pattern = "nodes\\w+")
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
  if (length(colnames(Nodes))>6){
    colnames(Nodes)=c('KW','category','node_weight','cls_cnt','btw_cnt','strength',colnames(Nodes)[seq(7,length(colnames(Nodes)))])  
  } else{
    colnames(Nodes)=c('KW','category','node_weight','cls_cnt','btw_cnt','strength')
  }
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
  edges$inv_edge_weight[edges$inv_edge_weight == 0] <- 1e-15 # vaules of 0 throw an error in Knet
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
    if (all(distances(net,V(net),V(net))==0)){
      nodes$cls_cnt <- 0
      }else{
        cls=closeness(net,normalized = T)
      nodes$cls_cnt <- cls[match(nodes$KW, names(cls))] 
      }
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
    qscore <- mapply(score,nodes$KW,nodes$btw_cnt,nodes$cls_cnt,nodes$strength)
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
  return(list('nodes'=Nodes,'edges'=edges))
}

t=lapply(files,SANTAsub)

# for (i in names(t)){
#   #print(i)
#   df=t[[i]]$nodes
#   #print(colnames(df))
#   write.table(df,file = i,row.names = F,col.names = T,quote = F,sep = '\t')
#   #print(paste(i,'.tsv',sep = '',collapse = ''))
# }
#print("SANTA: NODE RANKINGS ASSIGNED AND STORED WITHIN THE NODES TABLE")
#print("##########################################################################")
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
widnet=function(links,nodes,source,target,nodeid,edge_weight,col_weight,group,nodesize,click=T){
  links$V1=apply(as.matrix(links[,c(source)]),1,function(v)(
    as.numeric(rownames(nodes[nodes[,c(nodeid),drop=T]==v,]))-1))
  links$V2=apply(as.matrix(links[,c(target)]),1,function(v)(
    as.numeric(rownames(nodes[nodes[,c(nodeid),drop=T]==v,]))-1))
  rbPal <- colorRampPalette(c('black','red'))
  cuts=length(table(links[,c(col_weight),drop=T]))
  if (cuts>1){
    #cuts=length(table(links[,c(col_weight),drop=T]))
    links$col=rbPal(cuts)[as.numeric(cut(links[,c(col_weight),drop=T],breaks = cuts))] # coloured edges
  }
  else{
    links$col='#000000' #rbPal(1)[rep(1,nrow(links))] # all black
  }
  MyClickScript <- ifelse(isFALSE(click),NULL,
    sprintf('alert(d.name + " statistics: node size based on %s of " + d.nodesize + "; if varying, edge weights and colors represent the size of supporting literature for a given co-occurrence");',nodesize)
    )
  #
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
                 zoom = T,
                 Nodesize = nodesize,linkDistance = JS("function(d){return 2-Math.exp(1-d.value)}"),
                 radiusCalculation = "Math.pow(d.nodesize,1/2)+2",
                 opacity = 0.7, legend = TRUE, opacityNoHover = 0.3,
                 colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"),
                 linkColour=links$col,clickAction = MyClickScript)
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
  edge_weight='edge_weight',
  col_weight = 'R1C1',
  group='category',
  nodesize = 'qscore')
  saveWidget(p, file=paste(wd,tag,"_interactive_Gene-MeSH_Network.html",sep=''))
}
###
NetDist=function(edges){
  print("generate graph")
  g.tm <- graph_from_data_frame(edges,directed = FALSE) # conf as edge attribute
  #CC=components(g.tm)
  # use only 1 completely connected graph
  #CompNodes = which(CC$membership == 1)
  #g.tm = induced_subgraph(g.tm,CompNodes)
  #
  print(c("graph size:",gsize(g.tm)))
  ### DISTANCE MATRICES FOR KNET
  print("compute distance metrics...")
  D <- DistGraph(g.tm,edge.attr = "inv_edge_weight", dist.method = "shortest.paths", verbose=T) 
  B <- BinGraph(D,verbose=T,nsteps = 1000)
  return(list('graph'=g.tm,'B'=B, 'nodes'=V(g.tm)))
}
#
NetS=list("NEIGHnet"=NetDist(edges))
#
print("get cliques...")
cliqueset <- max_cliques(NetS$NEIGHnet$graph,min = 3,subset = sort(V(NetS$NEIGHnet$graph)$name))
print("...done")
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
  #
  set.seed(2202)
  result_knet <- Knet(g.tm, nperm=500, vertex.attr=pws, edge.attr="weight", verbose=T,B = B)
  #result_knet <- Compactness(g.tm, nperm=300, vertex.attr=pws, edge.attr="weight", verbose=T)
  if (length(pws) == 1){
    result_knet=list(result_knet)
    names(result_knet)=sprintf("%s",pws[1])
  }
  rm(g.tm,B,pos = environment())
  gc(verbose = F)
  return(result_knet)
}
#
for (data in names(NetS)){
  print("estimate available memory...")
  gc(verbose=F)
  free=as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE)) # available memory at that time
  use=0.25*free
  print(sprintf("set to use around 25%% of the available RAM, i.e. %g GB", round(use/1E6)))  
  net=NetS[data][[1]]
  g.tm=net$graph
  B=net$B
  # object.size() returns bytes
  # needs a 2x moltiplication factor because R parallel sucks!
  required=2*as.numeric(object.size(B)+object.size(g.tm)*500)/1E3 # returns KB
  print(sprintf("given B and 500 Knet permutations,each core will need %f GB",required/1E6))
  # set how many cores will then be used
  mcores=max(min(floor(use/required),ncores),1)
  print(sprintf("Then %g cores will be used",mcores))
  print("initialize parallelization")
  chunks=get.chunks(lpaths.in.graph,min(length(lpaths.in.graph),mcores))
  cl <- makeCluster(length(chunks),type = 'PSOCK')
  clusterEvalQ(cl, library(SANTA))
  clusterExport(cl, c("g.tm","B","Kstat"))
  print(c("Knet applied to",data))
  result_knet=parLapply(cl, chunks, Kstat)
  stopCluster(cl)
  print(c("Knet completed",data))
  print("Prepare results table...")
  result_knet=do.call("c", result_knet)
  names(result_knet)=str_replace(names(result_knet),"[0-9]+.","")
  pvals=sapply(result_knet, function(res){return(res$pval)})
  pvals=sort(pvals, decreasing = FALSE)
  enrichment<- data.frame(names(pvals), unname(pvals), p.adjust(unname(pvals), method = 'holm'), row.names = NULL,stringsAsFactors = F)
  colnames(enrichment)=c('cliqueID','pKnet','pAdjust')
  rm(result_knet)
  gc(verbose=F)
  print("...done")
}  
#
### select significant cliques ###
hit.cliques=lpaths.in.graph[enrichment[enrichment$pAdjust<=0.05,]$cliqueID]
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
#plot(g,layout=layout_nicely,vertex.size=4,vertex.label=NA)
g=induced_subgraph(g,names(hot.KW))
#
#plot(g,layout=layout_nicely,vertex.size=4,vertex.label=NA)
# EDGE WEIGHT HAS TO BE RE-SET!
g <- set.edge.attribute(g, "weight", index=E(g), 1-E(g)$inv_edge_weight)
#
### CLIQUES MIGHT BREAK IN COMMUNITY DETECTION! NEW CLIQUE-CLIQUE GRAPH BASED ON JACCARD (DISCONTINUED) ###
jaccard <- function(a, b){
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
#
# x=combn(hit.cliques,2,simplify=F,function(l){
#   return(data.frame(names(l)[1],names(l)[2],-log10(1-jaccard(l[[1]],l[[2]]))))
# })
# using geometric mean 
"geometric.mean" <- function(x,na.rm=TRUE){exp(mean(log(x),na.rm=na.rm))}
# using indicator function 
print("calculate intersections of cliques...")
cl <- makeCluster(ncores,type = 'PSOCK')
clusterExport(cl, c("g"))
### "x" consists of a collapsed distance measure between members of two cliques. Used to build a clique-clique graph ###
com <- vector(mode = "list", length = choose(length(hit.cliques),2))
com <-combn(hit.cliques,2,simplify=F)
x=vector(mode = "list", length = length(com))
x=parLapply(cl,com,function(l) {
  inter <- intersect(l[[1]],l[[2]])
  vec <- unique(l[[1]],l[[2]]) 
  Ivec <- 0+(vec %in% inter)#(l[[1]] %in% inter)+(l[[2]] %in% inter) #-length(inter) # indicator function
  df=data.frame(names(l)[1],names(l)[2],mean(c(igraph::E(g)[l[[1]] %--% l[[2]]]$weight,Ivec)),stringsAsFactors=F)
  colnames(df)=c('c1','c2','weight')
  return(df)
})
#
stopCluster(cl)
rm(com)
gc(verbose=F)
print("...done")
### immediately remove edges with 0 weight ###
x=x[sapply(x,function(df)(df$weight!=0))]
#
n=length(x)
dt <- data.table(c1=rep(NA,n), c2=rep(NA,n),weight=rep(0,n))
dt=rbindlist(x)
#colnames(x)=c('c1','c2','weight')
#wecdf=Vectorize(ecdf(x$weight))
#x$weight=-log10(1-x$weight)
#
### CLIQUE-CLIQUE GRAPH ###
print("generate cliques network...")
GC=graph_from_data_frame(dt,directed = F)
print("clean environment...")
rm(x,dt,pos=environment())
gc(verbose=F)
### PRUNE EDGES WITH 0 WEIGHT ###
#print("generate cliques network...")
GC=delete.edges(GC, which(E(GC)$weight < quantile(E(GC)$weight,0.5)))
#CC=cluster_louvain(GC) 
CC=cluster_leiden(GC,objective_function = 'modularity',beta=0,n_iterations=1000)# ,weights = NULL,resolution_parameter = 1.5,beta = 0.001,n_iterations = 1000)
#
#CC=cluster_leiden(g,objective_function = 'modularity',n_iterations = 1000)
#
if (length(communities(CC)) < 2){
  quit(save="no", status=13)
} else{
  print("check passed: at least 2 communities")
}
#
print(sprintf("the graph's modularity is %f",modularity(GC,weight=E(GC)$weight,membership = CC$membership)))
png(paste(tag,"_clique-clique_communities_net.png",sep=''),width = 12, height = 12,units = 'in',bg='grey95',res=200)
plot(GC, vertex.color=rainbow(max(CC$membership), alpha=0.6)[CC$membership],layout=layout.graphopt,edge.width=E(GC)$weight,vertex.label.cex=0.65)+title("communities in cliques-only graph")
dev.off()
### Propagate memberships to individual entities ### 
clust_named=communities(CC)
### sugiyama ###
layex <-  layout_with_sugiyama(GC, layers=apply(sapply(clust_named,
                                                       function(x) V(GC)$name %in% as.character(x)),
                                                1, which))
#
origvert <- c(rep(TRUE, vcount(GC)), rep(FALSE, nrow(layex$layout.dummy)))
realedge <- as_edgelist(layex$extd_graph)[,2] <= vcount(GC)
png(paste(tag,"_clique-clique_communities_Sugiyama.png",sep=''),width = 16, height = 12,units = 'in',bg='grey95',res=200)
plot(layex$extd_graph,
     edge.arrow.size=.1,
     edge.width=E(GC)$weight,vertex.label.cex=0.65,
     vertex.color=rainbow(max(CC$membership), alpha=0.7)[CC$membership],
     vertex.size=5,
     vertex.shape='circle',
     vertex.label=NA,
     sub='Use it to visualize connectivity between communities'
)+title("cliques-cliques Sugiyama graph")
dev.off()
rm(layex,pos=environment())
gc(verbose=F)
#### CLUSTER PROPAGATION #### 
clust_prop=lapply(clust_named,function(C)(unique(do.call('c',hit.cliques[C]))))
for (i in names(clust_prop)){
  clust_prop[[i]]=data.frame('entity'=clust_prop[[i]],'membership'=i,stringsAsFactors = F)
}
clust_prop=bind_rows(clust_prop)
clust_prop$membership=as.numeric(clust_prop$membership)#+1
# MULTIMAPPING WILL HAVE THEIR OWN MEMBERSHIP
print("resolve multi-mapping nodes...")
multi=as.character(unique(clust_prop[duplicated(clust_prop$entity),]$entity))
sizes=table(clust_prop$membership)
tabs=table(sizes)
print(sprintf("there are %g ties among community sizes", length(tabs[tabs>1])))
#print(tabs)
#print(multi)
### iterative assignment ###
for (i in multi){
  set.seed(2202)
  clust_prop[clust_prop$entity == i,]$membership=as.numeric(sample(
                                                            names(sizes[sizes==min(sizes[names(sizes) %in% clust_prop[clust_prop$entity==i,]$membership])])
                                                            ,1))
  #
  sizes=table(clust_prop$membership)
  #print(sizes)
}
#clust_prop[clust_prop$entity %in% multi,]$membership=sapply(clust_prop[clust_prop$entity %in% multi,]$entity,
#                                                            function(e){
cnodes=clust_prop
cnodes$nodesize=1
cnodes=cnodes[!(duplicated(cnodes)),]                                                              
### SAVE nodes-to-cliquecommunity table ###
write.table(clust_prop,paste(tag,"_nodes_in_clique-clique_communities_net.tsv",sep=''), sep="\t", row.names = F)
### RE-SUBSET co-occ net
g=NetS$NEIGHnet$graph
knodes_named=setNames(nodes$qscore,nodes$KW)
g=set_vertex_attr(g,"qscore", value = knodes_named[match(V(g)$name,names(knodes_named))])
g=induced_subgraph(g,cnodes$entity)
clust_prop=setNames(as.numeric(clust_prop$membership),clust_prop$entity)
### plot clique memberships propagated on entities
png(paste(tag,"_gene-MeSH_CliqueNetwork_by_clique_communities.png",sep=''),width = 12, height = 12,units = 'in',bg='grey95',res=200)
plot(g, vertex.color=rainbow(max(clust_prop), alpha=0.6)[unname(clust_prop[V(g)$name])],layout=layout.graphopt,edge.width=10*E(g)$weight,vertex.label.cex=0.65,sub='ambiguous nodes are assigned to the smallest communities')+title("Gene-MeSH, cliques-only network coloured by clique-based communities")
dev.off()
###
#### SAVE IT ALSO AS INTERACTIVE FORMAT #### 
print("generate interactive network of cliques...")
### convert igraph ###
V(g)$group=unname(clust_prop[V(g)$name])
g_toD3=igraph_to_networkD3(g, group = V(g)$group)
### transfer edge attributes ###
g_toD3$links$edge_weight=apply(g_toD3$links[,c(1,2)],1,function(v){
  v=as.character(g_toD3$nodes[v+1,]$name)
  return(E(g)[v[1] %--% v[2]]$edge_weight)
  })
#
g_toD3$links$R1C1=apply(g_toD3$links[,c(1,2)],1,function(v){
  v=as.character(g_toD3$nodes[v+1,]$name)
  return(E(g)[v[1] %--% v[2]]$R1C1)
  })
# keep nodesize to 1
g_toD3$nodes$nodesize=knodes_named[match(g_toD3$nodes$name,names(knodes_named))]
#
rbPal <- colorRampPalette(c('black','red'))
cuts=length(table(g_toD3$links$R1C1))
if (cuts>1){
  #cuts=length(table(links[,c(col_weight),drop=T]))
  g_toD3$links$col=rbPal(cuts)[as.numeric(cut(g_toD3$links$R1C1,breaks = cuts))] # coloured edges
}else{
  g_toD3$links$col='#000000'}
#
MyClickScript <- sprintf('alert(d.name + " statistics: node size based on %s of " + d.nodesize + "; if varying, edge weights and colors represent the size of supporting literature for a given co-occurrence");','nodesize')
#
p=forceNetwork(
  Links=g_toD3$links,
  Nodes=g_toD3$nodes,
  Source='source',
  Target='target',
  NodeID ='name',
  Value='edge_weight',
  linkColour = g_toD3$links$col,
  Group='group',
  Nodesize = 'nodesize',
  zoom=T,
  opacity = 0.7, legend = TRUE, opacityNoHover = 0.3,
  radiusCalculation = "Math.pow(d.nodesize,1/3)+2",
  colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"),clickAction=MyClickScript)
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
print("initialize parallelization...")
cl <- makeCluster(ncores,type = 'PSOCK')
clusterEvalQ(cl, library(igraph))
clusterExport(cl, c("knodes_named","max_attempts","combs","g"))
#
print("CONQUER!")
#
link=combn(okCC,2,simplify = F,function(l){
  # l is a pair of communities
  # le are the vertices that share an edge with each others' community members
  le=attr(E(g)[l[[1]] %--% l[[2]]],'vnames')
  if (length(le) > 0){
    lv=str_split(le,'\\|')
    clusterExport(cl, c("l"),envir=environment())
    print(sprintf("evaluating %i keywords-seeds linking communities %s...",length(lv),paste(names(l),collapse = ' and ')))
    # look at the neighbourhood of each seed
    CC_queues=parLapply(cl,X = lv,function(v){
      #print(v)
      nlv=V(g)$name[neighbors(g,v)]
      nlv=setdiff(nlv,v)
      # such neighbourhood must only contain members of the two communities
      nlv=nlv[nlv %in% c(l[[1]],l[[2]])]
      if (length(nlv)>combs-2){
        #print(nlv)
        queues=combn(nlv,combs-2,simplify=F,function(s)(c(v,unlist(s))))
        #scorele=E(g)[v[1] %--% v[2]]$weight*knodes_named[names(knodes_named) %in% v]
        #scorenle=lapply(nlv,function(n)(setNames(sum(E(g)[c(v,nlv) %--% n]$weight)*knodes_named[names(knodes_named)==n],n)))
        #rankings=combn(scorenle,combs-2,simplify=F,function(s)(c(scorele,unlist(s,use.names = T))))
        rankings=lapply(queues, function(q)(unlist(lapply(q,function(n)(setNames(sum(E(g)[q %--% n]$weight)*knodes_named[names(knodes_named)==n],n))))))
        print(rankings)
        #return(attr(V(g)[v],"Knode"))
        return(rankings)  
      }
      else{
        return(list())
      }      
    })
    CC_queues=sort(do.call('c',lapply(unlist(CC_queues,recursive = F),function(query)(setNames(sum(query),paste(sort(names(query)),sep='_',collapse = '_'))))),decreasing = T)
    return(CC_queues[!(duplicated(CC_queues))][1:max_attempts])  
  }
  else{
    return(setNames(rep(0,max_attempts),rep(NA,max_attempts)))
  }
})
#
print("CONQUERED!")
stopCluster(cl)
#
queries=do.call(rbind.data.frame, lapply(link,names))
colnames(queries)=seq(1,ncol(queries))
#
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}
#
queries=delete.na(queries,1)
#
write.table(queries,paste(wd,tag,"_ordered_queries.tsv",sep=''), sep="\t", row.names = F)
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