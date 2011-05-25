#!/usr/local/bin/Rscript

require(illuminaMousev2BeadID.db)

args <- commandArgs(trailingOnly=TRUE)
symbol <- args[1]

annot <- mget(ls(illuminaMousev2BeadIDSYMBOL), illuminaMousev2BeadIDSYMBOL, ifnotfound = NA)

res <- annot[grep(symbol, annot, ignore.case=T)]

if(length(res)>0){
res
}else{stop("No matches found")}
