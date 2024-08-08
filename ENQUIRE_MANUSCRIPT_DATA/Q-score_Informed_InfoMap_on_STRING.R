### cluster analysis on Q-score informed STRING network ###
library(igraph)
setwd('~')
qscoredgraph=read_graph("Ferroptosis_and_Immune_System_informed_STRING_PhysNet_with_QScores.graphml",format = 'graphml')
V(qscoredgraph)$qweight=Vectorize(ecdf(V(qscoredgraph)$weight))(V(qscoredgraph)$weight)
V(qscoredgraph)$probweight=(V(qscoredgraph)$qweight-median(V(qscoredgraph)$qweight))/(1-median(V(qscoredgraph)$qweight))
set.seed(2202)
clusts=cluster_infomap(qscoredgraph,e.weights = NULL,v.weights = V(qscoredgraph)$qweight,nb.trials = 1E3)
commies=communities(clusts)
### compare to unweighted infomap ### 
set.seed(2202) 
clustsunw=cluster_infomap(qscoredgraph,e.weights = NULL,v.weights = NULL,nb.trials = 1E3)
### Find communities that differ from unweighted infomap and with non-zero grand total ###
unwcommies=communities(clustsunw)
notinunw=parallel::mclapply(mc.cores = 16,commies,
                            FUN = function(i)(all(sapply(unwcommies,function(j)(length(setdiff(i, j)) > 0 || length(setdiff(i, j)) > 0)))))
newcommies=commies[unlist(notinunw)]
### Extract only communities with strictly positive grand total weight ###
newcommies=newcommies[which(sapply(newcommies,function(v)(sum(V(qscoredgraph)$weight[V(qscoredgraph)$name %in% v])))>0)]
### order by mean qscore
newcommies=newcommies[order(sapply(newcommies,function(v)(mean(V(qscoredgraph)$weight[V(qscoredgraph)$name %in% v]))),decreasing = T)]
commieonlynet=induced_subgraph(qscoredgraph,unique(unlist(newcommies[1:3])))
set.seed(2202)
l <- qgraph::qgraph.layout.fruchtermanreingold(get.edgelist(commieonlynet,names=FALSE),vcount=vcount(commieonlynet),niter=1e5,
                                               area=1*(vcount(commieonlynet)^2),repulse.rad=(vcount(commieonlynet)^1.8))
V(commieonlynet)$label.family="sans"
V(commieonlynet)$frame.color='gray'
E(commieonlynet)$color='white'
V(commieonlynet)$label.dist=.5
V(commieonlynet)$label.degree=-pi/2
V(commieonlynet)$label.color='black'
V(commieonlynet)$label.cex=.6
par(mar=c(0,0,0,0))
plot(commieonlynet,
     mark.groups = newcommies[1:3],
     vertex.size=ifelse(V(commieonlynet)$weight==0,0,3),
     vertex.label=ifelse(V(commieonlynet)$weight==0,NA,V(commieonlynet)$name),
     layout=l,
     vertex.color=ifelse(V(commieonlynet)$weight==0,'white','magenta3'))
#### ORA #### 

print("perform ORA to test overlap between Reactome pathways and the observed gene communities")

gs2ora=as.list(newcommies)
gs2ora=lapply(gs2ora,function(v)(v[v!='']))

# check gene sets with the same genes in #

gscollapse=sapply(gs2ora,function(v)(paste(sort(v),collapse = '_',sep='_')))

gsdups=sapply(unique(gscollapse),function(v)(which(unname(gscollapse)==v)))

gsdups=gsdups[sapply(gsdups,length)>1]

if (length(gsdups)>0){
  for (i in gsdups){
    gs2merge=unique(unlist(gs2ora[i]))
    gs2mergename=paste(names(gs2ora[i]),sep='_AND_',collapse = '_AND_')
    gs2ora[[gs2mergename]]=gs2merge
  }
  gs2ora[unname(unlist(gsdups))]=NULL
}

print(sprintf("A total %i unique genes are arranged into %i distinct communities",
            length(unique(unlist(gs2ora))),length(gs2ora)))


# # LOAD REACTOME GENE SETS ###
# 
# load(opt$NetS2)

# subset pathways according to the universe ()

paths2uni=pathways[names(pathways) %in% unique(unlist(gs2ora))]

okpaths=table(unname(paths2uni))
okpaths=names(okpaths)[okpaths>=2]

paths2uni=paths2uni[unname(paths2uni) %in% okpaths]

paths2uni=tapply(names(paths2uni),factor(unname(paths2uni)),function(v)(v))

print(sprintf("a total of %i pathways match the universe and with a resulting size bigger than 1",length(paths2uni)))

# also for Reactome Pathways we have to merge identical sets #

print("Check if any such universe-restricted pathways are identical and merge those...")

gscollapse=sapply(paths2uni,function(v)(paste(sort(v),collapse = '_',sep='_')))

gsdups=sapply(unique(gscollapse),function(v)(which(unname(gscollapse)==v)))

gsdups=gsdups[sapply(gsdups,length)>1]

if (length(gsdups)>0){
  for (i in gsdups){
    gs2merge=unique(unlist(paths2uni[i]))
    gs2mergename=paste(names(paths2uni[i]),sep='_AND_',collapse = '_AND_')
    paths2uni[[gs2mergename]]=gs2merge
  }
  paths2uni[unname(unlist(gsdups))]=NULL
}

print(sprintf("a total of %i Fisher's exact tests per observed gene set will be performed",length(paths2uni)))

print("ORA!")

print("note: null overlaps between observed gene sets and Reactome pathways are discarded without testing")

universe=unique(unlist(paths2uni))
gs2ora=lapply(gs2ora,function(v)(v[v %in% universe]))

gs2oratests=lapply(gs2ora,function(gs){
  gsc=setdiff(universe,gs)
  fishes=lapply(paths2uni,function(p){
    gsover=intersect(gs,p)
    R1C1=length(gsover)
    if (R1C1>0){
      R1C2=length(setdiff(gs,p))
      R2C1=length(intersect(gsc,p))
      R2C2=length(setdiff(gsc,p))
      ft=fisher.test(matrix(c(R1C1,R1C2,R2C1,R2C2),nrow = 2,ncol = 2,byrow = T),alternative = 'g')
      ft$overlap=paste(gsover,collapse = '_',sep='_')
      return(ft)
    }else{
      return(NULL)
    }
  })
  return(fishes[sapply(fishes,length)>0])
})

gs2oraresdf=data.frame('CommunityID'=unlist(mapply(function(x,y)(rep(x,y)),names(gs2oratests),sapply(gs2oratests,length),SIMPLIFY=F)),
                       'pathway'=unlist(sapply(gs2oratests,names,simplify=F)),
                       'ORA.P.Value'=unlist(sapply(gs2oratests,function(v)(sapply(v,'[[','p.value')),simplify=F)),
                       'ORA.Adj.P.Value'=p.adjust(unlist(sapply(gs2oratests,function(v)(sapply(v,'[[','p.value')),simplify=F)),'BH'),
                       'Community.Pathway.Overlap'=unlist(sapply(gs2oratests,function(v)(sapply(v,'[[','overlap')),simplify=F))
)

gs2oraresdf=gs2oraresdf[order(gs2oraresdf$ORA.P.Value),]
gs2oraresdf=gs2oraresdf[gs2oraresdf$ORA.Adj.P.Value<=0.01,]
sumscommie=lapply(newcommies,function(n)(mean(V(qscoredgraph)$weight[V(qscoredgraph)$name %in% n])))
gs2oraresdf$CommunityMeanQScore=unlist(sumscommie)[gs2oraresdf$CommunityID]
gs2oraresdf$CommunitySize=sapply(newcommies[gs2oraresdf$CommunityID],length)
gs2oraresdf$CommunityGenes=sapply(newcommies[gs2oraresdf$CommunityID],paste,sep='_',collapse='_')

write.table(gs2oraresdf,file = "Ferroptosis_and_Immune_System_Qscore-Informed_Infomap_Communities_and_ORA.tsv",sep='\t',quote=F,row.names = F)
