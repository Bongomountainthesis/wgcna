#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly=TRUE)

#optionally cutoff membership at a given pvalue
p.cut <-  args[1]

cat('plotting modules for p value cut-off of: ', p.cut,"\n")

load("results/E.RData")
load("results/WCGNA_part1.RData")
load("results/results.RData")

if(!file.exists("results/modules")) dir.create("results/modules")
for(i in unique(moduleLabels)){
  
  these <- moduleLabels==i
  nm <- unique(moduleColors[these])
  choose.cols <- c("ID","EnsemblID","symbol","chromosome_name","start_position","end_position", "strand", "description", paste("ME",nm, sep=""), paste("ME",nm,".p",sep=""))
  this.data <- results[these,choose.cols]

  #write the whole module out, sorted by p.value
  this.data<-this.data[order(this.data[,paste('ME',nm,'.p', sep="")]),]
  write.csv(this.data, paste("results/modules/", "module_",i,"_",nm,".csv",sep=""))
  
  if(!is.null(p.cut)){
    this.data <- this.data[this.data[,paste("ME",nm,".p",sep="")]<=p.cut,]
  }

  if(nrow(this.data)>0){
    this.xpn <- t(E[,this.data[,"ID"]])
    
    file=paste("results/modules/","module_",i,"_",nm , ".ps", sep="")
    
    postscript(file=file, paper="special", width=6, height=6)
    plot(this.xpn[1,], ylim=range(this.xpn), xlab="Sample", ylab="relative log2 intensity", type="l")

    if(nrow(this.xpn)>1){
      for(j in 2:nrow(this.xpn)){
        lines(this.xpn[j,],col=j, type="l")
      }
    }
    dev.off()  
  }
}
