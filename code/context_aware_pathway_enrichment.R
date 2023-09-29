#
if (isFALSE(exists("libloc"))) {
  libloc=.libPaths()[1]  
}
#
options(Ncpus = 4L)
print(c("library in use:",libloc))
### install missing libraries ###
print("check for missing, required packages...")
###### install libraries #######
if (!("BiocManager" %in% installed.packages(lib.loc=libloc)[,"Package"])){
  install.packages("BiocManager",repos="https://cloud.r-project.org/",lib=libloc,destdir = libloc)
}
list.of.packages <- c("optparse","data.table","RColorBrewer","openxlsx","igraph","stringr","dplyr","SANTA","parallel")
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
####### parse arguments ######
# create parser object

library(optparse)

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

option_list <- list( 
  make_option(c("-w", "--directory"), action="store", default=getwd(),
              help="Working directory (default to current working directory)",dest="wd",
              type="character",metavar="path"),
  make_option(c("-o", "--outdirectory"), action="store", default=getwd(),
              help="Output directory (default to current working directory)",dest="od",
              type="character",metavar="path"),
  make_option(c("-n","--netpathdata"),action="store", type="character", default='input/ENQUIRE-KNet_STRING_RefNet_Reactome_Paths.RData.gz', 
              help="Path to 'ENQUIRE-KNet_STRING_RefNet_Reactome_Paths.RData.gz' (required).\n\t\tIf the current working directory is not the 'ENQUIRE' folder, the default path will throw an error.",
              metavar="path",dest="NetS2"),
  make_option(c("-e", "--edgetable"),action="store", type="character", default=NULL, 
              help="Path to an ENQUIRE-generated, gene-gene edge table file (required)",
              metavar="path",dest="test_net"),
  make_option(c("-c", "--cores"),action="store", type="numeric", default=4, 
              help="max number of cores used (PSOCK parallelization) (default: 4), >1 recommended",
              metavar="parameter",dest="ncores"),
  make_option(c("-t", "--tag"),action="store", type="character", default="ENQUIRE", 
              help="tag prefix (default to 'ENQUIRE')",
              metavar="tag",dest="tag"),
  make_option(c("-s", "--setsize"),action="store", type="numeric", default=100, 
              help="maximum Reactome pathway size (default: 100, minimum 3)",
              metavar="parameter",dest="s"))

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults 
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$test_net)){
  print("ERROR: missing inputs!")
  print("Do 'Rscript context_aware_pathway_enrichment.R -h' for printing help")
  stop()
}
print("parsed options:")
print(str(opt))
print("#####")
###

library(SANTA)
library(stringr)
library(dplyr)
library(parallel)
library(igraph)
library(data.table)

####### preload ####
print("load precomputed data...")
load(file.path(opt$NetS2))
NetS2$STRINGnet$graph=upgrade_graph(NetS2$STRINGnet$graph)
g.ref=NetS2$STRINGnet$graph
####### 
print("#############################################")
tm.edges=read.delim(file = opt$test_net, stringsAsFactors = FALSE, header = T, row.names = NULL)
#
g.tm <- graph_from_data_frame(tm.edges,directed = FALSE) # conf as edge attribute
#
if (!is.null(E(g.tm)$inv_edge_weight)){
  g.tm=set_edge_attr(g.tm,'weight',index=E(g.tm),E(g.tm)$inv_edge_weight)
}else if(!is.null(E(g.tm)$ztPois.cdf)){
  g.tm=set_edge_attr(g.tm,'weight',index=E(g.tm),1-E(g.tm)$ztPois.cdf)
}else{
  print("WARNING: no valid edge weight found - will set edge weight to 1")
  g.tm=set_edge_attr(g.tm,'weight',index=E(g.tm),1)
}
test_nodes=intersect(V(g.ref)$name,V(g.tm)$name)
###
print(sprintf("found %i protein-coding, ENQUIRE-derived genes to map",length(test_nodes)))
###
degs=degree(g.tm)
tm.nodes=V(g.tm)$name
### for every pair of textmined v-v paths we compute 1/2^|d| ###
print("generate distance matrices...")
#
g2=g.ref
###
d.in=distances(g.ref,V(g.ref)$name,test_nodes,weights = NA,algorithm = c("dijkstra"))
d.tm=distances(g.tm,test_nodes,test_nodes,algorithm = c("dijkstra"),weights =E(g.tm)$weight)#,algorithm = c("dijkstra"))
# remove Inf distances
d.in[d.in==Inf]=max(d.in[d.in!=Inf])+1
d.tm[d.tm==Inf]=max(d.tm[d.tm!=Inf])+1
#
all.n=lapply(test_nodes,function(v)(list('START'=v,
                                         'N'=intersect(neighbors(g.tm,v)$name,test_nodes),
                                         'Nw'=d.tm[v,intersect(neighbors(g.tm,v)$name,test_nodes)],#setNames(E(g.tm)[v %--% intersect(neighbors(g.tm,v)$name,test_nodes)]$R1C1,intersect(neighbors(g.tm,v)$name,test_nodes)),
                                         'Nc'=intersect(setdiff(V(g.tm)$name,c(v,neighbors(g.tm,v)$name)),test_nodes),
                                         'Ncw'=d.tm[v,intersect(setdiff(V(g.tm)$name,c(v,neighbors(g.tm,v)$name)),test_nodes)])))
names(all.n)=test_nodes
#
Degs=degree(g.ref)
nsize=length(V(g.ref))
print("apply network propagation score Q on neighborhood...")
ncores=opt$ncores
cl <- makeCluster(ncores,type = 'PSOCK')
clusterExport(cl, c("d.in","d.tm","tm.nodes","test_nodes","g.ref","all.n","Degs","nsize"))
clusterEvalQ(cl,library(igraph))
####
v.ref5=parLapplyLB(cl,V(g.ref)$name,function(v){
  #
  print(v)
  #
  D=Degs[names(Degs)==v]
  #
  hm=(unlist(lapply(all.n[setdiff(names(all.n),v)],function(p){
    #print(p[['START']])
    ### COMPONENT 2.1: HOW MANY "CONSERVED NEIGHBOURS" FROM TXT TO STRING ### (arbitrary positive value)
    NfitN=sum(unlist(lapply(p[['N']][p[['N']] %in% colnames(d.in)],function(n){
      Nfit=as.numeric(d.in[v,p[['START']]]+d.in[v,n]<=2)*(exp(-p[['Nw']][names(p[['Nw']])==n])) #/length(p[['N']])
      ### COMPONENT 2.2: HOW LIKELY IS THIS VALUE, CONSIDERING THE NODE DEGREE AND SIZE OF TEXTMINED NET ### (CDF)
      #pN=ifelse(Nfit>0,dbinom(Nfit,D,pTxt),0)
      #print(Nfit)
      return(Nfit)
    })))
    NcfitN=sum(unlist(lapply(p[['Nc']][p[['Nc']] %in% colnames(d.in)],function(n){
      Nfit=as.numeric(d.in[v,p[['START']]]+d.in[v,n]<=2)*(exp(-p[['Ncw']][names(p[['Ncw']])==n])) #/length(p[['N']])
      ### COMPONENT 2.2: HOW LIKELY IS THIS VALUE, CONSIDERING THE NODE DEGREE AND SIZE OF TEXTMINED NET ### (CDF)
      #pN=ifelse(Nfit>0,dbinom(Nfit,D,pTxt),0)
      #print(Nfit)
      return(Nfit)
    })))
    return(NfitN+NcfitN)
  })))
  ##
  hmm=sum(hm)/D #mean(hm)
  S=hmm
  #
  if (is.nan(S)){
    return(setNames(0,v))
  }
  else{
    return(setNames(S,v))
  }
})
####
stopCluster(cl)
##
#### SOME PLOTTING ####
V(g2)$weight=0
#V(g2)[V(g2)$name %in% test_nodes]$weight=1
V(g2)$weight=unlist(v.ref5[match(names(unlist(v.ref5)),V(g2)$name)])
#### APPLY RULE OF THREE TO NEVER-WITNESSED CONNECTIONS IN PERMUTATION TESTS (10k permutations performed) #### 
E(g2)$weight=ifelse(E(g2)$weight==0,1E-8,E(g2)$weight)
#### NETWORK PRUNING (NECESSARY TO LOWER COMPUTATIONAL DEMANDS)
Isolated = which(degree(g2)==0)
g2 = delete.vertices(g2, which(degree(g2)==0))
### remove multiple edges and self loops
g2=delete.edges(g2,E(g2)[which_loop(g2)])
### QUANTILE NORMALIZATION OF vertex weights
qvref=Vectorize(ecdf(V(g2)$weight))
V(g2)$qweight=qvref(V(g2)$weight)
#V(g2)$qweight=V(g2)$krank
par(mfrow=c(1,1))
png(file.path(opt$od,paste0(opt$tag,"_Qscore_VS_degree.png")),width = 12, height = 12,units = 'in',bg='grey95',res=200)
plot(V(g2)$weight~degree(g2),xlab='degree',ylab='score')+title(sprintf("Correlation for %s scores: %f",opt$tag,cor(V(g2)$weight,degree(g2))))
dev.off()
#### APPLY KNET ####
### COMPUTE DISTANCE MATRIX ###
print("prepare matrices for Knet...")
#
#### COMPUTE PATHWAYS TO TEST ####ENQUIRE-KNet_STRING_RefNet_Reactome_Paths.RData.gz
tm.nodes=V(g.tm)$name
tm.in.pathways <- pathways#[names(pathways) %in% tm.nodes] # which pathways contain ANY textmined gene? 
okp=table(unname(tm.in.pathways))
###
### SIZE CONSTRAIN ###
okp=names(okp[okp > 2 & okp <= opt$s])
###
### FILTER: HIT SIZE > 2 AND AT LEAST 10% REPRESENTATIVENESS ###
#pw.ids = as.character(ptab[(ptab$inn > 2 & ptab$inn/ptab$all >= 0.1),]$path) #unique(unname(pathways[pathways %in% okp])) # which pathways contain AT LEAST 2 textmined genes? 
#pathways.in.graph=pathways[pathways %in% pw.ids] # then get ALL genes in these assessed pathways.
pw.ids = unique(unname(pathways[pathways %in% okp])) # which pathways contain AT LEAST 2 textmined genes? 
pathways.in.graph=pathways[pathways %in% pw.ids & names(pathways) %in% V(g2)$name]
### COMPUTE OVERLAP BETWEEN GENE SETS AND ACTUAL NODES IN GRAPH
over=unlist(lapply(unique(pathways.in.graph), function(pw)(setNames(length(intersect(names(pathways.in.graph[pathways.in.graph==pw]),names(pathways[pathways==pw])))/length(names(pathways[pathways==pw])),pw)))) # then get ALL genes in these assessed pathways
### "OPTIMAL STOP QUANTILE"
thr=quantile(over,exp(-1))
over=over[over>thr]
pathways.in.graph=pathways.in.graph[pathways.in.graph %in% names(over)]
### SET BOUNDARIES FOR A PATHWAY SET SIZE (KNET RESTRICTIONS) ###
ptab=table(unname(pathways.in.graph))
ptab=names(ptab[ptab > 2 & ptab < 100])
pw.ids = unique(unname(pathways.in.graph[pathways.in.graph %in% ptab]))
pathways.in.graph=pathways.in.graph[pathways.in.graph %in% pw.ids]
lpaths.in.graph=sapply(unique(pw.ids), function(pw){
  genes=names(pathways.in.graph[pathways.in.graph==pw])
  #genes=intersect(genes,tm.nodes)
  #genes=intersect(genes,V(NetS2$STRINGnet$graph)$name)
  names(genes)=rep(pw,length(genes))
  return(genes)})
#
lpaths.in.graph=lpaths.in.graph[sapply(lpaths.in.graph, function(v)(length(v)>2))]
#
get.chunks=function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
#
chunks=get.chunks(lpaths.in.graph,ncores)
#
### custom knet (more efficient)
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
###
Kstat=function(chunk){
  #
  pws=names(chunk)
  badp=c()
  exp.nodes=V(g.tm)$name
  for (pw in pws) {
    if (sum(as.numeric(exp.nodes %in% unname(chunk[pw][[1]]))) >= 3){
      pheno=as.numeric(exp.nodes %in% unname(chunk[pw][[1]]))*V(g.tm)$weight
      if (sum(pheno)==0){
        badp=c(badp,pw)
      }
      else{
        g.tm <- set_vertex_attr(
          g.tm, name=pw,
          value=pheno)
      }
    }
    else{
      badp=c(badp,pw)
    }
    # g.tm <- set_vertex_attr(
    #       g.tm, name=pw,
    #       value=V(g.tm)$weight)
  }
  #
  pws=setdiff(pws,badp)
  #
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
###
for (data in names(NetS2)){
  print("initialize parallelization")
  cl <- makeCluster(ncores,type = 'PSOCK')
  clusterEvalQ(cl, library(SANTA))
  net=NetS2[data][[1]]
  g.tm=net$graph
  B=net$B
  clusterExport(cl, c("g.tm","B","Kstat","customKnet"))
  clusterEvalQ(cl,environment(customKnet) <- asNamespace('SANTA'))
  clusterEvalQ(cl,assignInNamespace("Knet", customKnet, ns = "SANTA"))  
  print(c("Knet applied to",data))
  result_knet=parLapply(cl,chunks, Kstat)
  stopCluster(cl)
  gc()
  print(c("Knet completed",data))
  print("Prepare results table...")
  result_knet=do.call("c", result_knet)
  names(result_knet)=str_replace(names(result_knet),"[0-9]+.","")
  #
  #result_knet <- Knet(g.tm, nperm=100, vertex.attr=pw.ids, edge.attr="conf", verbose=T,B = B,parallel = 32)
  pvals=sapply(result_knet, function(res){return(res$pval)})
  pvals=sort(pvals, decreasing = FALSE)
  png(file.path(opt$od,paste0(opt$tag,"_KNet_pvalue_distribution.png")),width = 12, height = 12,units = 'in',bg='grey95',res=200)
  hist(pvals)
  dev.off()
  ### APPLY BH BASED ON PRIORS DISTRIBUTION INSTEAD OF U(0,1) #### 
  qvals=p.adjust(pvals,method = 'holm')
  critq=length(qvals[qvals<=0.05])
  print(c("CRITICAL Q-VALUE AT POSITION",critq))
  ###
  print("save output...")
  #fraction <- vector(mode = "character", length = length(names(pvals)))
  percentage <- vector(mode = "numeric", length = length(names(pvals)))
  pathway_size <- vector(mode = "numeric", length = length(names(pvals)))
  overlap <- vector(mode = "numeric", length = length(names(pvals)))
  i <- 1
  for (pw in names(pvals)) {
    pathway_size[i] <- sum(unname(pathways)==pw)
    overlap[i] <- length(lpaths.in.graph[pw][[1]])
    #fraction[i] <- paste("[",paste(in.network, total, sep="/"),"]", sep = "")
    percentage[i] <- round((overlap[i]/pathway_size[i])*100,2) 
    i <- i+1
  }
  #
  cats=pathcat[match(names(pvals),pathcat)]
  #cats=cats[!duplicated(cats)]
  #
  enrichment<- data.frame(names(cats),names(pvals), unname(pvals), qvals, p.adjust(unname(pvals), method = 'holm'),sapply(result_knet[names(pvals)],'[[','AUK.obs'), pathway_size,overlap,percentage, row.names = NULL,stringsAsFactors = F)
  colnames(enrichment) <- c("category","pathway", "p_KNet","adj_p_modBH", "adj_p_Holm","AUK.obs","pathway_size","overlap","%")
  write.table(enrichment,file.path(opt$od,paste0(opt$tag,"_KNet_pathway_enrichment.tsv")),row.names = F,quote = F, sep = '\t')
  #results[sprintf("%s",data)]=list(enrichment)
}
###
