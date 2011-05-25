#!/usr/local/bin/Rscript

require("WGCNA")

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)

probe <- as.character(args[1])
cut.off <- as.numeric(args[2])

probe
cut.off



load("results/E.RData")
load("results/WCGNA_part1.RData")
load("results/results.RData")
load("data/TOM.RData")


these <- TOM[probe,]
other.probes <- rownames(TOM)[which(these>=cut.off)]
other.probes<-other.probes[-1*which(other.probes == probe)]
these.results <- results[other.probes,]
dim(these.results)
write.csv(these.results, file=paste(probe,"_nearbytom.csv", sep="")) 

##########get TOM of gene of interest for cytoscape

#remove NAs

these.results <- these.results[which(these.results[,"symbol"]!="NA"),]

ids <- these.results[,"ID"]

nodeNames <- these.results[,"symbol"]

probeTOM <- TOM[probe,ids]

probeTOM <- as.data.frame(probeTOM)
                                        #Export TOM to cytoscape
cat("Exporting for Cytoscape\n")

cyt = exportNetworkToCytoscape(
  probeTOM,
  nodeNames = nodeNames,
  edgeFile = "Probe-CytoscapeInput-edges-FullTOM.txt",
  nodeFile = "Probe-CytoscapeInput-nodes-FullTOM.txt",
  weighted = TRUE,
  threshold = 0.5,
  nodeAttr = moduleColors
  )




