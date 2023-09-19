### fuzzy c means ### 
system("export DOWNLOAD_STATIC_LIBV8=1")
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("At least three arguments must be supplied (main dir, working dir, tag)", call.=FALSE)
}else{
  #print(args)
  sd=args[1]
  wd=args[2]
  tag=args[3]
  if (length(args)==5){
    edge_tab=args[4]
    node_tab=args[5]  
  }
}
#
if (isFALSE(exists("libloc"))) {
  libloc=.libPaths()[1]  
}
#
options(Ncpus = 8L)
print(c("library in use:",libloc))
### install missing libraries ###
print("check for missing, required packages...")
#install.packages(c(),repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc)
###
# factoextra and randomcolor gridBase excluded
if (!("BiocManager" %in% installed.packages(lib.loc=libloc)[,"Package"])){
  install.packages("BiocManager",repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc)
}
list.of.packages <- c("RColorBrewer","igraph","BiocManager","stringr","ggplot2","dplyr","qgraph","dynamicTreeCut","ppclust","svglite","gridBase")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages(lib.loc=libloc)[,"Package"])]
if (length(new.packages)>=1){  
  BiocManager::install(new.packages,lib=libloc,destdir = libloc) # FAU mirror
}
### handle failed installations ###
failed.packages <- new.packages[!(new.packages %in% installed.packages(lib.loc=libloc)[,"Package"])]
print("FAILED PACKAGES:")
print(failed.packages)
print("##########")
if (length(failed.packages)>=1){
  ### pre-install ad hoc dependencies ###
  packageurl='https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.2.1.tar.gz'
  install.packages(packageurl, repos=NULL, type='source',lib=libloc,destdir = libloc) # FAU mirror  
  packageurl='https://cran.r-project.org/src/contrib/Archive/lattice/lattice_0.20-40.tar.gz'
  install.packages(packageurl, repos=NULL, type='source',lib=libloc,destdir = libloc) # FAU mirror  
  packageurl='https://cran.r-project.org/src/contrib/Archive/statmod/statmod_1.4.35.tar.gz'
  install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
  packageurl="https://cran.r-project.org/src/contrib/Archive/latticeExtra/latticeExtra_0.6-29.tar.gz"
  install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
  packageurl="https://cran.r-project.org/src/contrib/Archive/lme4/lme4_1.1-32.tar.gz"
  install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
  packageurl="https://cran.r-project.org/src/contrib/Archive/acepack/acepack_1.4.0.tar.gz"
  install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
  packageurl <- "https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-75.tar.gz"
  install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc) 
  packageurl = "https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.3-1.tar.gz"
  install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc) 
  packageurl='https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-52.tar.gz'
  install.packages(packageurl, repos=NULL, type='source',lib=libloc,destdir = libloc) # FAU mirror  
  packageurl <- "https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.4-4.tar.gz" #install.packages(packageurl, repos=NULL, type="source")
  install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
  packageurl="https://cran.r-project.org/src/contrib/Archive/survival/survival_3.2-7.tar.gz"
  install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
  ### proceed with libraries of interest ###
  install.packages(failed.packages,repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc) # FAU mirror
}

###
# packageurl='https://cran.r-project.org/src/contrib/Archive/statmod/statmod_1.4.35.tar.gz'
# install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
# packageurl="https://cran.r-project.org/src/contrib/Archive/lme4/lme4_1.1-25.tar.gz"
# install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.4-4.tar.gz" #install.packages(packageurl, repos=NULL, type="source")
# install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.4-4.tar.gz" #install.packages(packageurl, repos=NULL, type="source")
# install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
# install.packages("car",repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc)
# install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
# install.packages("qgraph",repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc)
# packageurl = "https://cran.r-project.org/src/contrib/Archive/FactoMineR/FactoMineR_2.4.tar.gz"
# install.packages(packageurl, repos=NULL, type='source',lib=libloc,destdir = libloc) # FAU mirror  
# install.packages("factoextra",repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc)
# packageurl='https://cran.r-project.org/src/contrib/Archive/curl/curl_4.3.1.tar.gz'
# install.packages(packageurl, repos=NULL, type='source',lib=libloc,destdir = libloc) # FAU mirror  
# packageurl='https://cran.r-project.org/src/contrib/Archive/V8/V8_3.3.1.tar.gz'
# install.packages(packageurl, repos=NULL, type='source',lib=libloc,destdir = libloc) # FAU mirror  
# install.packages(failed.packages,repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc)
# 1) car recipe
# if (any(c("pbkrtest","car") %in% failed.packages)){
#   packageurl <- "https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.4-4.tar.gz" #install.packages(packageurl, repos=NULL, type="source")
#   install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
#   install.packages("car",repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc)
# }
# # 2) qgraph recipe
# if (any(c("foreign","Hmisc","qgraph") %in% failed.packages)){
#   packageurl <- "https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-75.tar.gz"
#   install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
#   packageurl = "https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.3-1.tar.gz"
#   install.packages(packageurl, repos=NULL, type="source",lib=libloc,destdir = libloc)
#   install.packages("qgraph",repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc)
# }
# # 3) factoextra recipe
# if (any(c("FactoMineR","factoextra") %in% failed.packages)){ 
#   packageurl = "https://cran.r-project.org/src/contrib/Archive/FactoMineR/FactoMineR_2.4.tar.gz"
#   install.packages(packageurl, repos=NULL, type='source',lib=libloc,destdir = libloc) # FAU mirror  
#   install.packages("factoextra",repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc)
# } 

# # 4) ggplot2 recipe 
# if (any(c('ggplot2') %in% failed.packages)){
#   packageurl='https://cran.r-project.org/src/contrib/Archive/lattice/lattice_0.20-40.tar.gz'
#   install.packages(packageurl, repos=NULL, type='source',lib=libloc,destdir = libloc) # FAU mirror  
#   packageurl='https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-52.tar.gz'
#   install.packages(packageurl, repos=NULL, type='source',lib=libloc,destdir = libloc) # FAU mirror  
#   install.packages(failed.packages,repos="https://ftp.fau.de/cran/",lib=libloc,destdir = libloc)
# }


### FUZZY CLUSTERING ANALYSIS ####
library(stringr)
library(ggplot2)
library(dplyr)
library(parallel)
library(igraph)
library(data.table)
library(ppclust)
library(gridBase)
library(grid)
library(RColorBrewer)
#
###
setwd(wd)
#
if (all(c(exists('node_tab'),exists('edge_tab')))){
  test_net <- read.delim(paste(wd,"/",edge_tab,collapse = '',sep=''),stringsAsFactors = F)
  test_net <- read.delim(paste(wd,"/",node_tab,collapse = '',sep=''),stringsAsFactors = F)
}else{
  test_net <- read.delim(paste(wd,"/",tag,"_Complete_edges_table_subgraph.tsv",collapse='',sep=''),stringsAsFactors = F)
  test_nodes <- read.delim(paste(wd,"/",tag,"_Complete_nodes_table_subgraph.tsv",collapse='',sep=''),stringsAsFactors = F)
}
net=graph_from_data_frame(test_net,directed = F,vertices = test_nodes)
### CENTROIDS ARE THE MESH TERMS ### 
centroids=test_nodes[,1,drop=T][grepl("MESH",test_nodes[,2,drop=T])]
## use invlogweighted similarity ##
sims=matrix(similarity(net,method='invlogweighted'),#method='invlogweighted'),
            ncol=vcount(net),
            nrow=vcount(net),
            dimnames = list('colnames'=V(net)$name,'rownames'=V(net)$name))
###
diag(sims)=Vectorize(function(n)(sum(rep(1,length(neighbors(net,n)))/log(degree(net)[neighbors(net,n)$name]))))(V(net)$name)
### compute ECDF to standardize data ### 
#sims[is.numeric(sims)]=Vectorize(ecdf(sims))(sims[is.numeric(sims)])
### make it so that it's a similarity measure (smaller=closer)
#sims[is.numeric(sims)]=1-sims[is.numeric(sims)]
### centmatrix = meshXmesh ###

centmat=sims[centroids,]
testmat=sims[-match(centroids,rownames(sims)),]
### use ward clustering to simplify genes sets - at least 4 kws each 
wardtree=hclust(dist(centmat,method = 'manhattan'),method = 'ward.D2')
wardtree$height <- round(wardtree$height, 6) 
clustmesh=setNames(dynamicTreeCut::cutreeDynamic(wardtree,minClusterSize = 4,deepSplit = 1,distM = as.matrix(dist(centmat,method = 'manhattan'))),wardtree$labels)
print("WARD AND DYNAMIC TREE CUT RESULTS:")
print(table(clustmesh))
### extract condensed centers by specifying memberships from ward/treecut ###
### also known as sharding ###
shards=tapply(names(clustmesh),factor(unname(clustmesh)),function(rows)(apply(centmat[rows,],2,mean)))
shards=do.call('rbind',shards)
###
u=inaparc::imembrand(nrow(sims), k=nrow(shards))$u
v=inaparc::kmpp(shards,k=nrow(shards))$v
#v=fcm(centmat,centers = length(centroids),fixcent = T,m=1.5,iter.max=2)
###
res.fcm <- fcm(sims,stand = F,centers=v,m=2,memberships = u,fixcent = F,nstart=5,iter.max = 1000,con.val=1e-9)
genes=match(rownames(testmat),rownames(sims))
#
#### use PCA then custom plot ####
pcadat=prcomp(rbind(v,sims),center = T)
####
### use defuzzyfied clusters to draw ellipses 
clusts=setNames(unname(res.fcm$cluster),rownames(sims)[as.numeric(names(res.fcm$cluster))])
pcashard=tapply(names(clusts),factor(unname(clusts)),function(v)(as.data.frame(matrix(pcadat$x[v,c('PC1','PC2')],
                                                                                      ncol = 2,
                                                                                      nrow = length(v),
                                                                                      dimnames = list(rownames=v,colnames=c('PC1','PC2'))))))
pcashardf=bind_rows(pcashard)

pcashardf$Cluster=clusts[match(rownames(pcashardf),names(clusts))] #unlist(sapply(names(pcashard),function(n)(rep(n,nrow(pcashard[[n]])))))
pcashard=pcashardf
#pcashard=cbind(data.frame('Cluster'=unlist(sapply(names(pcashard),function(n)(rep(n,nrow(pcashard[[n]]))))),row.names = names(pcashard)),
#               bind_rows(pcashard))
### genes set of at least 2 ###
okclust=levels(factor(pcashard$Cluster))[tapply(rownames(pcashard),factor(pcashard$Cluster),function(g)(sum(as.numeric(g %in% rownames(testmat)))))>=2]
pcashard=pcashard[pcashard$Cluster %in% okclust,]
pcashard$Cluster=factor(as.numeric(pcashard$Cluster),levels = sort(as.numeric(unique(pcashard$Cluster))))
##### implement cluster labels ####
### use the top meshs by Wscore to annotate ellipses ###
pcacent=tapply(names(clusts),factor(unname(clusts)),function(v)(as.data.frame(pcadat$x[v,c('PC1','PC2')])))
pcacent=cbind(data.frame('Cluster'=unlist(sapply(names(pcacent),function(n)(rep(n,nrow(pcacent[[n]])))))),
              bind_rows(pcacent))
### remove genes ###
pcacent=pcacent[!(rownames(pcacent) %in% rownames(testmat)),]
pcacent$sorting=test_nodes[match(rownames(pcacent),test_nodes$KW),ncol(test_nodes)-2]
# top meshs by membership define descriptions #
tops=res.fcm$u[!(rownames(res.fcm$u) %in% rownames(testmat)),]
tops=apply(tops,2,function(cc)(rownames(tops)[order(cc,decreasing = T)[1:5]]))
tops=lapply(colnames(tops),function(x)(paste(str_trunc(str_remove_all(tops[,x,drop=T],"\\/.+$"),width = 15), sep = '',collapse=', ')))
labs=setNames(unlist(tops),str_match(colnames(res.fcm$u),"[0-9]+"))
### non-abbreviated version ###
topslong=res.fcm$u[!(rownames(res.fcm$u) %in% rownames(testmat)),]
topslong=apply(topslong,2,function(cc)(rownames(topslong)[order(cc,decreasing = T)[1:5]]))
topslong=lapply(colnames(topslong),function(x)(paste(topslong[,x,drop=T], sep = '',collapse=', ')))
labslong=setNames(unlist(topslong),str_match(colnames(res.fcm$u),"[0-9]+"))
#x=pcacent %>% mutate(row=rownames(.)) %>% group_by(Cluster) %>% slice_max(order_by = sorting, n = 3,with_ties = F)
#labs=tapply(x$row,factor(x$Cluster),function(x)(paste(str_trunc(str_remove_all(x,"\\/.+$"),width = 20), sep = '',collapse=', ')))
pcacent$label=labs[match(pcacent$Cluster,names(labs))]
#pcacentcoord=as.data.frame(pcadat$x[grepl("Cl",rownames(pcadat$x)),c('PC1','PC2')])
#pcacentcoord$label=labs[as.numeric(str_match(rownames(pcacentcoord),"[0-9]+"))]
pcashard$lab=labs[match(pcashard$Cluster,names(labs))]
pcashard$lablong=labslong[match(pcashard$Cluster,names(labslong))]
#### save also complete mesh-cluster associations ####
fullmesh2cl=data.frame(cl=unname(res.fcm$cluster),en=rownames(res.fcm$u))
fullmesh2cl=fullmesh2cl[!(fullmesh2cl$en %in% rownames(testmat)), ]
fullmeshs=tapply(fullmesh2cl$en,factor(fullmesh2cl$cl),function(v)(paste(v,sep=', ',collapse = ', ')))
pcashard$MeSHFullSet=fullmeshs[pcashard$Cluster]
clust2name=setNames(pcashard$Cluster,pcashard$lab)
clust2name=clust2name[!(duplicated(names(clust2name)))]
###
### general purpose color vector ###
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#clustcols=list(setNames(randomcoloR::distinctColorPalette(k=ncol(res.fcm$u)),as.numeric(str_match(colnames(res.fcm$u),"[0-9]+"))))
clustcols=list(setNames(rep(col_vector,ceiling(ncol(res.fcm$u)/length(col_vector)))[1:ncol(res.fcm$u)],as.numeric(str_match(colnames(res.fcm$u),"[0-9]+"))))
clustcols[[1]]=setNames(ifelse(names(clustcols[[1]]) %in% okclust,clustcols[[1]],'#FFFFFF'),names(clustcols[[1]]))
clustcolleg=setNames(clustcols[[1]],pcashard[match(as.character(names(clustcols[[1]])),pcashard$Cluster),]$lab)
clustcolleg=clustcolleg[!duplicated(clustcolleg)]
names(clustcolleg)=(ifelse(is.na(names(clustcolleg)),'Other MeSHs',names(clustcolleg)))
#clustcolleg[names(clustcolleg)!='Other MeSHs']=randomcoloR::distinctColorPalette(k=length(clustcolleg)-1,runTsne = T)
#clustcolleg[names(clustcolleg)!='Other MeSHs']=randomcoloR::distinctColorPalette(k=length(clustcolleg)-1,runTsne = T)
clustcolleg[names(clustcolleg)!='Other MeSHs']=rep(col_vector,ceiling((length(clustcolleg)-1)/length(col_vector)))[1:(length(clustcolleg)-1)]
clustcols[[1]][okclust]=clustcolleg[names(clust2name)[match(okclust,clust2name)]]
names(clustcolleg)=tools::toTitleCase(names(clustcolleg))
clustcols[[1]][is.na(clustcols[[1]])]="#FFFFFF"
###
#pcashard$Genes=rownames(pcashard)
### 
write.table(pcashard[intersect(rownames(pcashard),rownames(testmat)),],paste(wd,"/",tag,"_GeneSets_FCM.tsv",collapse='',sep=''),row.names = T,col.names = NA,quote = F,sep = '\t')
###
##### PLOTTING
### make ad hoc legend ###
leg=ggplot() + geom_point(data=data.frame(x=runif(length(clustcolleg)),y=runif(length(clustcolleg)),
                                          Description=factor(names(clustcolleg),
                                                             levels=c(sort(make.unique(names(clustcolleg))[!grepl("Other MeSHs",names(clustcolleg))]),
                                                                      names(clustcolleg)[grepl("Other MeSHs",names(clustcolleg))]))),
                          aes(x=x,y=y,fill=Description),
                          shape=21,size=4,stroke=1,color='black')+
  scale_fill_manual(values=clustcolleg)+
  xlim(-1,-1)+ylim(-1,-1)+theme_void(base_size = 12)+theme(text = element_text(size=12,face='bold'),
                                                           legend.text.align = 0,
                                                           legend.spacing.y = unit(0.005,'in'),
                                                           legend.spacing.x = unit(0.05,'in'),
                                                           legend.position = 'bottom',
                                                           legend.box.just = 'bottom')+
  guides(fill=guide_legend(ncol=2,byrow=TRUE))
leg
### PIE CHART REPRESENTATION ###
props=lapply(rownames(res.fcm$u),function(r)(setNames(res.fcm$u[r,,drop=T],colnames(res.fcm$u))))
names(props)=rownames(res.fcm$u)



### account for vectors of just 0s

if (all(sapply(props,function(v)(all(v<=0.05))))){
  print("WARNING: ALL PREDICTED MEMEBERSHIP DEGREES ARE < 5%")
  print("The clusters are highly fuzzy")
  print("proportions will be kept intact")
}else{
  print('simplify memberships based on a minimum of 5% memberships')
  props=lapply(props,function(v)(ifelse(v >= .05,v,0)))
  props[sapply(props,function(v)(all(v < 0.05)))]=lapply(props[sapply(props,function(v)(all(v < 0.05)))],function(v)(rep(1/length(v),length(v))))  
}

l <- qgraph::qgraph.layout.fruchtermanreingold(get.edgelist(net,names=FALSE),vcount=vcount(net),
                                               area=10*(vcount(net)^2),repulse.rad=(vcount(net)^2.7))
# Set plot layout
svglite::svglite(paste(wd,"/",tag,"_GeneSets_FCM_plot.svg",collapse='',sep=''),bg = 'transparent',fix_text_size = T,width = 16,height = 18)
layout(mat = matrix(c(2, 1, 1, 2), 
                    nrow = 2, 
                    ncol = 2),
       heights = c(1,10),    # Heights of the two rows
       widths = c(1, 10))     # Widths of the two columns
par(mar=c(0,0,0,0))
plot.igraph(net,layout=l,
            vertex.shape='pie',
            vertex.pie=props,
            vertex.size=ifelse(test_nodes[match(V(net)$name,test_nodes[,1]),2]=='MESH',0.5,8),
            vertex.label=ifelse(test_nodes[match(V(net)$name,test_nodes[,1]),2]=='MESH',NA,V(net)$name),
            vertex.label.cex=ifelse(test_nodes[match(V(net)$name,test_nodes[,1]),2]=='MESH',0,1.5),
            vertex.pie.color=clustcols,
            label.x=1,
            label.y=-1,
            vertex.label.font=2)
## the last one is the current plot
plot.new()              ## suggested by @Josh
vps <- baseViewports()
pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
#plot(leg)
# create an apporpriate viewport.  Modify the dimensions and coordinates as needed
vp.BottomRight <- plotViewport()#c(5,0,0,0))#margins =c(1.8,1,0,1))
###
# plot the ggplot using the print command
print(leg,vp = vp.BottomRight)      
###         
dev.off()

