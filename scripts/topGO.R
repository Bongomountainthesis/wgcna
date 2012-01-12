#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly=TRUE)

#optionally cutoff membership at a given pvalue
p.cut <-  args[1]

stringsAsFactors=FALSE

load("results/E.RData")
load("results/WCGNA_part1.RData")
load("results/results.RData")

library(topGO)

data <- read.csv(file = "module_2_blue.csv")

geneList <- data[,"X"]


