#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  cat("\nDESCRIPTION:\nGenerates a graph that contains nodes for each of the genes in any of the chose modules, whose correlation with the eigengene of that module has a p-value of less than <max p>. Edges are drawn between these nodes if the topological overlap between the genes they represent is greater than <min TOM>. Edges are weighted by topological overlap.
\nUSAGE: igraph_pretty.R <min TOM> <max p> <module> <module> ...
ARGS:
     <min TOM> - the minimum topological overlap for which to draw an edge
     <max p> - the maximum p.value of correlation to the module eigengene for which to draw a node (gene)
     <module> - The names (colours) of the modules you want to plot\n\n
")
  options( show.error.messages=FALSE)
  stop()
}
 
library(beadarray)
library(WGCNA)
library(igraph)
options(stringsAsFactors = FALSE);

load("E.RData")
load("WCGNA_part1.RData")
load("results.RData")
load("adjacency.RData")
load("TOM.RData")

min.tom <- as.numeric(args[1])
max.p <- as.numeric(args[2])

mod.col <- as.character(args[c(-1,-2)])

plot.module <- function(mod.col, min.tom, max.p, TOM, wgcna.res, layout=layout.spring){

  colnm <- paste("ME",mod.col, sep="")

  #which genes are in these modules?
  these.i <- rownames(wgcna.res)[which(wgcna.res[,"ModuleName"] %in% mod.col)]

  #which of these are correlated with their module eigengene with a sig of <= max.p?
  colnm.p <- paste(colnm, 'p',sep=".")
  ok.p<-unique(unlist(  
	           lapply(colnm.p, 
                          function(x){
                             cat(x,"\n")
                             rownames(wgcna.res)[which(wgcna.res[,x] <= max.p)]
                          }
                   )
              )
         )
  these.i <- intersect(these.i, ok.p)

  if(length(these.i) ==0) stop('No members of modules with eigengene correlations of < max.p')

  #get the TOM data for the genes we have selected
  this.clust <- TOM[these.i,these.i]

                                        #this was supposed to stop it drawing a load of unconnected points, but as we're filtering on
                                        #pvalue now, I think we can get away with it.
                                        #  valid.rows <- rownames(this.clust)[apply(this.clust, 1, function(x){any(x[x!=1]>=min.tom)})]
                                        #  this.clust <- this.clust[valid.rows,valid.rows]
  
                                        #these.i <- intersect(valid.tom.rows,these.i) 
  
  if(nrow(this.clust)>0){
    
    ids <- rownames(this.clust)
    gr <- graph.adjacency(this.clust, 
                          mode="undirected", 
                          weighted=TRUE, 
                          diag=FALSE, 
                          add.colnames=TRUE)
    
    V(gr)$color <- wgcna.res[ids,"ModuleName"]
    
                                        #remove edges that have < min.tom connectivity
    rm <- which(E(gr)$weight < min.tom)
    rm <- rm-1  # edges are 0-based
    gr <- delete.edges(gr, E(gr)[rm])
    
    nms <-  wgcna.res[ids,"symbol"]
    inds <- which(nms=="NA")
    nms[inds] <- ids[inds] 
    V(gr)$label <- nms
    inds<-inds-1
    #matt wants unknown probes removed
    gr<-delete.vertices(gr, inds)
    
# possible layouts.
#     layout.random
#     layout.circle
#     layout.sphere
#     layout.fruchterman.reingold
#     layout.kamada.kawai
#     layout.spring
#     layout.reingold.tilford
#     layout.fruchterman.reingold.grid
#     layout.lgl
#     layout.graphopt
#     layout.mds
#     layout.svd
#     layout.norm

     

    
    #title=paste("Module", mod.col, "\n", "corr p.val", max.p, "topo overlap", min.tom)
    title=""
   plot(gr, layout=layout, main=title)
    
  }else{stop("No connections above min.tom")}
  
}



e.genes <- unique(gsub("ME","",gsub("\\.p","",grep("ME", colnames(results), value=T))))

if(length(mod.col)==0) mod.col <- e.genes

f <- paste(paste("modules/moduleGraph", paste(mod.col, collapse="_"), sep="_"), '.ps', sep="")

postscript(file=f, paper="special", height=6, width=6,  fonts=c("serif", "Palatino"), horizontal=FALSE)
plot.module(mod.col=mod.col, min.tom=min.tom, max.p=max.p, TOM=TOM, wgcna.res=results)
dev.off()

cat("\nDone. Results in file",f,"\n\n")
