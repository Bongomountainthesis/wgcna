#!/usr/local/bin/Rscript

library(beadarray)
library("illuminaMousev2BeadID.db")

BSData <- get(load("data/BSData.quantile.RData"))

args <- commandArgs(trailingOnly=TRUE)
probe <- as.character(args[1])

E<- exprs(BSData)
 #2^range(E)

#2^(E[probe,])
range(E)

E[probe,]
