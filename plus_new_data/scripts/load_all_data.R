############ Load, QC and normalise ###############

library(beadarray)

dataFile = "data/All_Samples_Unnormalised_25012011.txt"

BSData <- readBeadSummaryData(dataFile=dataFile,
                              skip=0,
                                ProbeID="ProbeID",
                              columns = list(exprs = "AVG_Signal",
                                se.exprs="BEAD_STDERR",
                                NoBeads = "Avg_NBEADS",
                                Detection="Detection Pval"
                                             ),
                              )
 
#> colnames(exprs(BSData))
# [1] "MGE1_160810"  "MGE2_160810"  "MGE3_160810"  "MGE4_160810"  "MGE5_160810"
# [6] "MGE6_160810"  "MGE8_160810"  "MGE9_160810"  "MGE1_300610"  "MGE2_300610"
#[11] "MGE3_300610"  "MGE4_300610"  "MGE5_300610"  "MGE6_300610"  "MGE7_300610"
#[16] "MGE8_300610"  "P4_RNA1"      "2_CHx6_P1"    "P10_RNA4"     "4_CHx6_P1"
#[21] "12_CHx6_P1"   "1_CHx6_P1"    "11_CHx6_P1"   "P10_RNA3"     "E12_5_13_sub"
#[26] "E13_5_19_sub" "6_sub"        "7_sub"        "E12_5_4_sub"  "5_sub"
#[31] "P4_RNA4_sub"  "E13_5_23_sub" "E13_5_21"     "E12_5_20"     "E12_5_12"
#[36] "E13_5_20"     "10"           "E13_5_25"     "8"            "4"
#[41] "E12_5_13"     "E13_5_19"     "6"            "7"            "E12_5_4"
#[46] "5"            "P4_RNA4"      "E13_5_23"     "9_CHx6_P1"    "P21_RNA3"
#[51] "P10_RNA1"     "P4_RNA3"      "P21_RNA4"     "10_CHx6_P1"   "5_CHx6_P1"
#[56] "E13_5_22"     "7_CHx6_P1"    "P10_RNA2"     "6_CHx6_P1"    "3_CHx6_P1"
#[61] "P21_RNA2"     "P21_RNA1"     "8_CHx6_P1"    "P4_RNA2"      "E12_5_19"
#[66] "E12_5_2"      "2"            "E12_5_18"     "3"            "E12_5_11"
#[71] "E12_5_5"      "E13_5_6"      "E13_5_3"      "E14_5_7"      "E12_5_10"
#[76] "E13_5_5"      "9"            "11"           "1"            "E13_5_2"

##get relevant arrays plus _sub arrays for ones that failed before (5, 6 and 7)
BSData <- BSData[,c("MGE1_160810",  "MGE2_160810" , "MGE3_160810" , "MGE4_160810" , "MGE5_160810","MGE6_160810" , "MGE8_160810" , "MGE9_160810",  "MGE1_300610" , "MGE2_300610","MGE3_300610" , "MGE4_300610" , "MGE5_300610" , "MGE6_300610" , "MGE7_300610","MGE8_300610","E12_5_13_sub", "E12_5_4_sub", "E12_5_20", "E12_5_12", "E12_5_13", "E12_5_4", "E12_5_19", "E12_5_5", "E12_5_18", "E12_5_11", "E12_5_5", "E12_5_10", "E13_5_19_sub", "E13_5_23_sub", "E13_5_21","E13_5_20", "E13_5_25", "E13_5_19", "E13_5_23", "E13_5_22", "E13_5_6", "E13_5_3", "E14_5_7", "E13_5_5", "E13_5_2")]

#BSData<-BSData[,mice]
save(BSData, file="data/BSData.RData")


# some QC plots:

postscript(file="results/Boxplotprenorm.ps", horizontal=FALSE)
boxplot(as.data.frame(log2(exprs(BSData))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,2,2,2,4,4,4,4,6,6,6), main="Log2 Expression")
dev.off()

postscript(file="results/plotMAXYprenorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData), arrays = 1:36,  pch = 16)
dev.off()

postscript(file="results/plotDENSITYprenorm.ps", horizontal=FALSE)
E <- log2(exprs(BSData))
plot(density(E[,1]))
for(i in 1:11){
  lines(density(E[,i]),col=i)
}
dev.off()

#From the QC plots, clearly some of these haven't worked:
fail <- c("E12_5_13_sub", "E12_5_4_sub","E13_5_19_sub", "E13_5_23_sub", "E12_5_11")
BSData<-BSData[, !colnames(exprs(BSData)) %in% fail ]

# Quantile Normalise
BSData.quantile = normaliseIllumina(BSData, method = "quantile", transform = "log2")


postscript(file="results/Boxplotpostnorm.ps", horizontal=FALSE)
boxplot(as.data.frame((exprs(BSData.quantile))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,2,2,2,4,4,4,4,6,6,6), main="Log2 Expression")
dev.off()

postscript(file="results/plotMAXYpostnorm_for_E115.ps", horizontal=FALSE)
plotMAXY(exprs(BSData.quantile), arrays = 1:16, log = FALSE, pch = 16)
dev.off()

save(BSData.quantile, file="data/BSData.quantile.RData")

postscript(file="results/plotDENSITYpostnorm.ps", horizontal=FALSE)
E <- exprs(BSData.quantile)
plot(density(E[,1]))
for(i in 1:10){
  lines(density(E[,i]),col=i)
}
dev.off()


################ Limma #################

library(limma)
library(biomaRt)

# Bioconductor packages are based on Mark Dunning's ReMOAT reannotations.
# The BeadID package uses numeric ArrayAddressIDs as keys and is for data 
# that has been summarised from bead level - so this is the file you need 
# if you are using the beadarray library to read in the raw bead data.
library(illuminaMousev2BeadID.db) 

# This package has much the same annotation, again based on the ReMOAT process
# but with ProbeIDs (ILMN_*) as keys. This is for data that has already been 
# collapsed into a single value per probe. Which is what we have in this case.
#
# Note that crap probes which don't map to anything or cross map to loads of
# stuff have been removed before annotation.
library(illuminaMousev2.db) 

E <- exprs(BSData.quantile)


# We only have 3 timepoints, which isn't *really* enough to treat the
# data as a proper timecourse, so just run Limma on all three possible
# pairwise comparisons.
design<-matrix(0,nrow=(ncol(E)), ncol=3)
colnames(design) <- c("Control","GF","GF_TNF")
rownames(design) <- colnames(E)
design[1:4,1] <- 1
design[5:8,2] <- 1
design[9:11,3] <- 1

cont.matrix<-makeContrasts(ControlvGF=GF-Control, 
			   ControlvGF_TNF=GF_TNF-Control,
                           GFvGF_TNF=GF_TNF-GF, 
                           levels=design)


fit<-lmFit(E, design)
fit<-contrasts.fit(fit, cont.matrix)
ebFit<-eBayes(fit)

ids = rownames(E)

symbol <- mget(ids, illuminaMousev2BeadIDSYMBOL, ifnotfound = NA)
sum(sapply(symbol, length) > 1)
symbol <- as.character(symbol)
length(which(symbol=="NA"))

ensembl = mget(ids, illuminaMousev2BeadIDENSEMBL, ifnotfound = NA)
length(which(is.na(ensembl)))

sum(sapply(ensembl, length) > 1)
crosshyb <- which(( sapply(ensembl, length) ) > 1)
length(crosshyb)
ensembl[crosshyb] <- NA
ensembl <- as.character(ensembl)
ensembl[ensembl=="NA"] = NA
length(which(is.na(ensembl)))

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

ebFit$genes = anno


# We're really doing a 1-way ANOVA here, with a post-hoc test between all possible
# comparisons.  The Ftest will tell you which differences are significant, so these are the
# p-values and adjusted p-values you are interested in.
# The contrasts can then tell you where the significant difference *is*. 
f.test<- topTable(ebFit, number=nrow(E))
rownames(f.test)<-f.test$ID
f.test<-f.test[order(f.test[,"P.Value"]),]
write.csv(f.test,"results/f_test.csv",row.names=F)

#ControlvGF=GF-Control
#ControlvGF_TNF=GF_TNF-Control
#GFvGF_TNF=GF_TNF-GF
#P4vP10=P10-P4,   P10vP21=P21-P10,       P4vP21=P21-P4

GF_control <- topTable(ebFit, coef=1, adjust="BH", number=nrow(E))
rownames(GF_control) <- GF_control$ID
GF_control <- GF_control[order(abs(GF_control[,"logFC"]), decreasing=TRUE),]
GF_control <- GF_control[!duplicated(GF_control[,"symbol"]),]
GF_control <- GF_control[order(GF_control[,"adj.P.Val"],decreasing=FALSE),]
write.csv(GF_control,"results/GF_control.csv",row.names=F)

TNF_control<-topTable(ebFit, coef=2, adjust="BH", number=nrow(E))
rownames(TNF_control)<-TNF_control$ID
TNF_control <- TNF_control[order(abs(TNF_control[,"logFC"]), decreasing=TRUE),]
TNF_control <- TNF_control[!duplicated(TNF_control[,"symbol"]),]
TNF_control <- TNF_control[order(TNF_control[,"adj.P.Val"],decreasing=FALSE),]
write.csv(TNF_control,"results/TNF_control.csv",row.names=F)

TNF_GF<-topTable(ebFit, coef=3, adjust="BH", number=nrow(E))
rownames(TNF_GF)<-TNF_GF$ID
TNF_GF <- TNF_GF[order(abs(TNF_GF[,"logFC"]), decreasing=TRUE),]
TNF_GF <- TNF_GF[!duplicated(TNF_GF[,"symbol"]),]
TNF_GF <- TNF_GF[order(TNF_GF[,"adj.P.Val"],decreasing=FALSE),]
write.csv(TNF_GF,"results/TNF_GF.csv",row.names=F)












