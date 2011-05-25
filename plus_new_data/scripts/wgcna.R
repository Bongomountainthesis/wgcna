#!/usr/local/bin/Rscript

require("beadarray")
require("genefilter")
require("WGCNA")
require("illuminaMousev2BeadID.db")
require("biomaRt")
require("lme4")
require("lattice")
require("illuminaMousev2.db") 

cat("loading data\n")

                                        #load the data
options(stringsAsFactors = FALSE);
BSData <- get(load("data/BSData.quantile.RData"))
E <- exprs(BSData)

                                        #remove anything that isn't expressed at more than the mean
                                        #expression level in at least 4 arrays.
                                        #This is pretty arbitrary

#cut.off <- quantile(E,0.75)
#cut.off <- mean(E)
cut.off <- 8
xpd <- genefilter(E,filterfun(function(x){sum(x >= cut.off) >= 5}))
sum(xpd)
E <- E[xpd,]

                                        #and remove anything remaining that is invariant
varying <- genefilter(E,filterfun(function(x){sd(x)>= 0.25}))
sum(varying)
E <- E[varying,]

                                        #WGCNA expects samples as rows, genes as cols
E <- t(E) 
save(E, file="results/E.RData") 


                                        #Fetch the Ensembl IDs and Gene Symbol from the Bioconductor Annotation
cat("fetching annotation from library\n")
ids = rownames(exprs(BSData))
symbol <- mget(ids, illuminaMousev2BeadIDSYMBOL, ifnotfound = NA)
#just check they only have one each
symbol <- as.character(symbol)

                                        #Also get Ensembl ID from BioC annotation so we can map through biomaRt to any other annotation we require:
ensembl = mget(ids, illuminaMousev2BeadIDENSEMBL, ifnotfound = NA)

                                        #Some ensembl genes which have multiple probes hitting them. For the moment, I am just looking at them seperately.
                                        #There are also a few probes which cross-hybridise to multiple ensembl transcripts. I have thrown their annotation away for now.
crosshyb <- which(( sapply(ensembl, length) ) > 1) 
ensembl[crosshyb] <- NA 
ensembl <- as.character(ensembl)
ensembl[ensembl=="NA"] = NA

                                        #Fetch other annotation from biomaRt
cat("fetching annotation from Biomart\n")
ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

filters <- "ensembl_gene_id"
values <- ensembl[!is.na(ensembl)]
attributes <- c("ensembl_gene_id", "chromosome_name","start_position", "end_position", "strand", "description")
ens.anno <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)
rownames(ens.anno)<-ens.anno[,1]

anno <- data.frame(
                   ID = as.character(ids),
                   EnsemblID = ensembl,
                   symbol=symbol,
                   ens.anno[ensembl,],
                   stringsAsFactors=F
              )
rownames(anno) <- anno[,"ID"]
save(anno, file="data/anno.RData")

#calculate range of soft threshold powers
#powers <- c(c(1:10), seq(from = 12, to=20, by=2))
#sft <- pickSoftThreshold(E, powerVector = powers, verbose = 5)

#Plot the results

#sizeGrWindow(9, 5)
#par(mfrow = c(1,2));
#cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
#postscript(file="results/timecourse_sft_threshold.ps",horizontal=FALSE)
#plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#main = paste("Scale independence"));
#text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
#abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
#plot(sft$fitIndices[,1], sft$fitIndices[,5],
#xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#main = paste("Mean connectivity"))
#text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()

#####################################################################################################
# Actual Work Starts HERE!

                                        #soft power settings
cat("Calculating adjacencies\n")
#softPower = 7;
#adjacency = adjacency(E, power = softPower, type="signed");
#save(adjacency, file="data/adjacency.RData")
load("data/adjacency.RData")
                                        #Calculate the TOM from the adjacency and the dist for the clustering from that. 
cat("Calculating TOM\n")
#TOM = TOMsimilarity(adjacency, TOMType="signed");
load("data/TOM.RData")
dissTOM = 1-TOM
#colnames(TOM) <- rownames(TOM) <- colnames(dissTOM) <- rownames(dissTOM)<- colnames(E)
#save(TOM, file="data/TOM.RData")

                                        #And do the clustering with their fast hclust variant
cat("Clustering\n")
geneTree = flashClust(as.dist(dissTOM), method = "average");

                                        # Plot the resulting clustering tree (dendrogram)
postscript(file="results/Dendro1.ps",paper="special",width=12,height=9, horizontal=F)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04);
dev.off()
                                        #Define some modules using dynamic tree cut
minModuleSize = 70;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize );

                                        #Give clusters colour names instead of numbers:

dynamicColors = labels2colors(dynamicMods)

                                        #And plot the dendrogram with the colours underneath

postscript(file="results/Dendro2.ps",paper="special", width=8,height=6, horizontal=F)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

                                        #Calculate eigengenes and cluster them
cat("Clustering Eigengenes\n")
MEList = moduleEigengenes(E, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = flashClust(as.dist(MEDiss), method = "average");

                                        #Plot the result
cls <- gsub("ME","",colnames(MEs))
file <- paste("results/MEs", 1, ".ps", sep="")  
postscript(file=file, paper="special", width=10, height=6, horizontal=F)
par(mfrow=c(3,1))
for(i in 1:ncol(MEs)){
    plot(MEs[,i], type="l", col=cls[i], xaxt="n", yaxt="n")
    if(i%%3==0){
      if(i != ncol(MEs)){
        dev.off()
        file <- paste("results/MEs", i+1, ".ps", sep="")  
        postscript(file=file, paper="special", width=10, height=6, horizontal=F)
        par(mfrow=c(3,1))
      }
    }
}            
dev.off()  


postscript(file="results/EG_clust.ps", paper="special", width=7, height=6, horizontal=F)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres = 0.4
abline(h=MEDissThres, col = "red")
dev.off()
                                        #Merge similar clusters on the basis of their eigengene similariy
cat("Merging clusters\n")
merge = mergeCloseModules(E, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;

                                        #Compare the new, merged, modules with the old ones
postscript(file='results/DendroMerge.ps', paper="special", width=12, height=9, horizontal=F)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

cls <- gsub("ME","",colnames(mergedMEs))
file <- paste("results/mergedMEs", 1, ".ps", sep="")  
postscript(file=file, paper="special",width=10, height=6, horizontal=F)  
par(mfrow=c(3,1))
for(i in 1:ncol(mergedMEs)){
    plot(mergedMEs[,i], type="l", col=cls[i], xaxt="n", yaxt="n")
    if(i%%3==0){
      if(i != ncol(mergedMEs)){
        dev.off()
        file <- paste("results/mergedMEs", i+1, ".ps", sep="")  
        postscript(file=file, paper="special", width=10, height=6, horizontal=F)
        par(mfrow=c(3,1))
      }
    }
}            
dev.off()  

mergedMEDiss = 1-cor(mergedMEs);
mergedMETree = flashClust(as.dist(mergedMEDiss), method = "average");
postscript(file="results/merged_EG_clust.ps", paper="special", width=7, height=6, horizontal=F)
plot(mergedMETree, main = "Clustering of merged module eigengenes", xlab = "", sub = "")
dev.off()

                                        #Save the merged clusters
moduleColors = mergedColors
colorOrder = c("grey", standardColors());
moduleLabels = match(moduleColors, colorOrder)-1;
names(moduleColors) <- names(moduleLabels) <- colnames(E)
MEs = mergedMEs;

                                        # Save module data
cat("Exporting Modules\n")
save(MEs, moduleLabels, moduleColors, geneTree, file = "results/WCGNA_part1.RData")


                                        #Make a dataframe of genes, their module memberships and their annotation 
                                        #Pearson cor of genes with eigengenes as module membership 
geneModuleMembership = as.data.frame(cor(E, MEs, use = "p"));
                                        #and associated p.value of that correlation. 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(E)));
colnames(MMPvalue) <- paste(colnames(MMPvalue), '.p', sep="")

geneids <- rownames(geneModuleMembership)
results <- cbind(anno[geneids,],ModuleNumber=moduleLabels[geneids], ModuleName=moduleColors[geneids], geneModuleMembership, MMPvalue[geneids,])

write.csv(results, file="results/results.csv", row.names=F)
save(results, file="results/results.RData")

    
                                        #Export TOM to cytoscape
#cat("Exporting for Cytoscape\n")
#nodeNames <- colnames(TOM)
#altNodeNames <- anno[nodeNames,"symbol"]
#cyt = exportNetworkToCytoscape(
#  TOM,
# nodeNames = nodeNames,
#  altNodeNames = altNodeNames,
#  edgeFile = "CytoscapeInput-edges-FullTOM.txt",
#  nodeFile = "CytoscapeInput-nodes-FullTOM.txt",
#  weighted = TRUE,
#  threshold = 0.8,
#  nodeAttr = moduleColors
#  )

                                        #Export TOME to VisANT
#cat("Exporting for VisANT")
#vis = exportNetworkToVisANT(TOM,
#  file = "VisANTInput-FullTOM.txt",
#  weighted = TRUE,
#  threshold = 0,
#  probeToGene = data.frame(nodeNames, altNodeNames )
#)


cat("Done!\n")
