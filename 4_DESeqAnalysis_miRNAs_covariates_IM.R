# Martijs Jonker / Ioannis Moustakas
# Perform a DEG analysis on miRNA
# 08-02-2016
# R-3.2.1
# cn-03
# /zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1209-Paul_Lucassen/MAD1209-P001-brain_tissue/MAD1209-P001-E001_2014_FFPE_miRNASeq_svleeuw1


########### Selectie maken
## miRNA's die een count < 20 hebben
# 30-10-2015

rm(list=ls())
options(stringsAsFactors = FALSE)
X11.options(type="Xlib")

design <- read.delim("./ExperimentInfo/SampleInformation/design.txt",stringsAsFactors=FALSE)
#df <- read.delim("./Results/miRNA_analysis_sel/CountTable.raw.sel.txt",stringsAsFactors=FALSE)
df <- read.delim("./Results/miRNA_analysis_sel/CountTable.raw.sel.txt",stringsAsFactors=FALSE)

# log2 transform the table ans save it

logTransformed <- log2(countsOnly+1)
row.names(logTransformed) <- df[,1] 
write.table(logTransformed, "./Results/miRNA_analysis_sel/CountTable.logTrans.sel.txt", sep="\t", quote=FALSE)

geneIDs <- df[,c(1)]
CountTable <- data.matrix(df[,-c(1)])
length(which(rowSums(CountTable) == 0)) # [1] 0

sizeFactors.mad <- function (counts, locfunc = median){
    loggeomeans <- rowMeans(log(counts))
    apply(counts, 2, function(cnts) exp(locfunc((log(cnts) - 
        loggeomeans)[is.finite(loggeomeans)])))
}
sf <- sizeFactors.mad(CountTable)
#countdata delen door de sizefactors#
CountTable.norm <- CountTable
for(i in 1:ncol(CountTable.norm)){CountTable.norm[,i] <- CountTable.norm[,i]/sf[i]}

# log2 transform the table ans save it
logTransformed <- log2(CountTable.norm + 1)
row.names(logTransformed) <- df[,1] 
write.table(logTransformed, "./Results/miRNA_analysisIoannis/CountTable.logTransSFNormalized.sel.txt", sep="\t", quote=FALSE)

boxplot(data.frame(log2(as.matrix(CountTable)+1)),pch=".", main="raw")
X11()
boxplot(data.frame(log2(as.matrix(CountTable.norm)+1)),pch=".", main="Norm")

colour <- c("red", "blue")
for(i in 1:ncol(CountTable.norm)){
  if(i == 1){
    plot(density(log2(as.matrix(CountTable.norm[ ,i])+1)), ylim=c(0, 0.3), col=colour[design$Group[i]])
  }
  if(i > 1){
    lines(density(log2(as.matrix(CountTable.norm[ ,i])+1)), col=colour[design$Group[i]])
  }
}

####
#Perform a PCA analysis showing the effects of the covariates Age, AgeGroup, Gender, NBB
#Perform the following tests:
# Depression VS control  within males
# Depression VS control  within females
# Depression VS control  within Young
# Depression VS control  within Old

library(annotate)
library(DESeq)

rownames(CountTable) <- geneIDs
cds <- newCountDataSet(CountTable,as.factor(design$Group))

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))
cds <- estimateDispersions(cds)
plotDispEsts(cds)

#adjust here#
cds@phenoData@data$condition
res <- nbinomTest(cds, unique(design$Group)[2], unique(design$Group)[1])
# head(res)
  
completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "Depr.vs.Cont_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()
  
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))

## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "Depr.vs.Cont_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()
  
completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "Depr.vs.Cont_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)


################# $$$$$$$$$$$$$$$$$$$$$ ####################
##### Try gender as group and see if there are DEGs there

rownames(CountTable) <- geneIDs
cds <- newCountDataSet(CountTable,as.factor(design$Gender))

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))
cds <- estimateDispersions(cds)
plotDispEsts(cds)

#adjust here#
cds@phenoData@data$condition
res <- nbinomTest(cds, unique(design$Gender)[2], unique(design$Gender)[1])
# head(res)

completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "Gender_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()

minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))

## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "Gender_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()

completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "Gender_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)


################# $$$$$$$$$$$$$$$$$$$$$ ####################
##### Try NNB (time sample was taken) as group and see if there are DEGs there

### First we need to group the NNB values into two groups. Do:
yearsSampledSorted <- sort(sapply(strsplit(design$NBB, "-"), function(x) x[1]))
# "1990" "1993" "1995" "1995" "1996" "1996" "1997" "1999" "2001" "2001" "2002" 
# "2005" "2005" "2006" "2006" "2008" "2008" "2009" "2009" "2010"
median(as.numeric(yearsSampledSorted))
# 2001.5 
# best to group as: Old, if <=2002, new if >2002 

# build a OLD/NEW vector and append to design
oldNewVector <- sapply(strsplit(design$NBB, "-"), function(x) x[1])<=2002
oldNewVector <- gsub(T, "OLD", oldNewVector)
oldNewVector <- gsub(F, "NEW", oldNewVector)
design$oldNew <- oldNewVector


rownames(CountTable) <- geneIDs
cds <- newCountDataSet(CountTable,as.factor(design$oldNew))

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))
cds <- estimateDispersions(cds)
plotDispEsts(cds)

#adjust here#
cds@phenoData@data$condition
res <- nbinomTest(cds, unique(design$oldNew)[2], unique(design$oldNew)[1])
# head(res)

completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "SampleAge_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()

minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))

## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "SampleAge_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()

completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "SampleAge_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)


# completeTablePath = paste0("./Results/miRNA_analysisIoannis/DesignTable.txt")
# write.table(design, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)


################# $$$$$$$$$$$$$$$$$$$$$ ####################
##### Try ageGroup as group and see if there are DEGs there

rownames(CountTable) <- geneIDs
cds <- newCountDataSet(CountTable,as.factor(design$AgeGroup))

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))
cds <- estimateDispersions(cds)
plotDispEsts(cds)

#adjust here#
cds@phenoData@data$condition
res <- nbinomTest(cds, unique(design$AgeGroup)[2], unique(design$AgeGroup)[1])
# head(res)

completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "AgeGroup_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()

minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))

## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "AgeGroup_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()

completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "AgeGroup_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)


#Perform the following tests:
# Normalize data
# Make PCAs
# Depression VS control  within sample OLD
# Depression VS control  within sample NEW
# Depression VS control  within males
# Depression VS control  within females
# Depression VS control  within Young
# Depression VS control  within Old
# Save maplots, number DEGs with Pvalue < 0.0 5 and adjPvalue < 0.05#


############ $$$$$$$$$$$$$$ ############### 
designAllFactors <- read.delim("./Scripts/TestBasicQC/seqdesign.txt", stringsAsFactors=FALSE)

# Depression VS control within sample OLD

# slice the design and count table and keep Old samples only
oldSamplesDesign <- designAllFactors[designAllFactors$sampleAge=="OLD", ]
oldSampleNames <- oldSamplesDesign$SampleID
oldSamplesCountTable <- CountTable[,oldSampleNames]
rownames(oldSamplesCountTable) <- geneIDs

#redo the DEG analysis as before
cds <- newCountDataSet(oldSamplesCountTable, as.factor(oldSamplesDesign$sampleinfo))

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))
cds <- estimateDispersions(cds)
plotDispEsts(cds)

#adjust here#
cds@phenoData@data$condition
res <- nbinomTest(cds, unique(oldSamplesDesign$sampleinfo)[2], unique(oldSamplesDesign$sampleinfo)[1])
# head(res)

completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "oldSamples_Dep.vs.Cont_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()

minPval <- min(res$pval[!is.na(res$pval)])
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))


## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "oldSamples_Dep.vs.Cont_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()

completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "oldSamples_Dep.vs.Cont_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)

# minPval
minPval

# minPadj
minPadj

# geneCount
geneCount



############ $$$$$$$$$$$$$$ ############### 
designAllFactors <- read.delim("./Scripts/TestBasicQC/seqdesign.txt", stringsAsFactors=FALSE)

# Depression VS control within sample NEW

# slice the design and count table and keep Old samples only
oldSamplesDesign <- designAllFactors[designAllFactors$sampleAge=="NEW", ]
oldSampleNames <- oldSamplesDesign$SampleID
oldSamplesCountTable <- CountTable[,oldSampleNames]
rownames(oldSamplesCountTable) <- geneIDs

#redo the DEG analysis as before
cds <- newCountDataSet(oldSamplesCountTable, as.factor(oldSamplesDesign$sampleinfo))

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))
cds <- estimateDispersions(cds)
plotDispEsts(cds)

#adjust here#
cds@phenoData@data$condition
res <- nbinomTest(cds, unique(oldSamplesDesign$sampleinfo)[2], unique(oldSamplesDesign$sampleinfo)[1])
# head(res)

completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "newSamples_Dep.vs.Cont_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()

minPval <- min(res$pval[!is.na(res$pval)])
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))


## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "newSamples_Dep.vs.Cont_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()

completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "newSamples_Dep.vs.Cont_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)

# minPval
minPval

# minPadj
minPadj

# geneCount
geneCount

############## $$$$$$$$$$$$$$$ #################
# none of the above works (no DEGs). So use the PCA plot for sample age and separate the samples according to PC1. 
# Group one above 0 (the less degraded) and Group two, PC1 below 0, (the more degraded) 

# group one 
lessDegraded = c("S03", "S04", "S05", "S08", "S09", "S10", "S12", "S18", "S19", "S20")

# group two
moreDegraded = c("S01", "S02", "S06", "S07", "S11", "S13", "S14", "S15", "S16", "S17")


############ $$$$$$$$$$$$$$ ############### 
designAllFactors <- read.delim("./Scripts/TestBasicQC/seqdesign.txt", stringsAsFactors=FALSE)

# Depression VS control within sample NEW

# slice the design and count table and keep less Degraded
lessDegradedDesign <- designAllFactors[designAllFactors$SampleID %in% lessDegraded, ]
lessDegradedNames <- lessDegradedDesign$SampleID
lessDegradedCountTable <- CountTable[,lessDegradedNames]
rownames(lessDegradedCountTable) <- geneIDs

#redo the DEG analysis as before
cds <- newCountDataSet(lessDegradedCountTable, as.factor(lessDegradedDesign$sampleinfo))

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))
cds <- estimateDispersions(cds)
plotDispEsts(cds)

#adjust here#
cds@phenoData@data$condition
res <- nbinomTest(cds, unique(lessDegradedDesign$sampleinfo)[2], unique(lessDegradedDesign$sampleinfo)[1])
# head(res)

completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "lessDegraded_Dep.vs.Cont_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()

minPval <- min(res$pval[!is.na(res$pval)])
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))


## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "lessDegraded_Dep.vs.Cont_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()

completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "lessDegraded_Dep.vs.Cont_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)
lessDegradedRes <- res

# minPval
minPval

# minPadj
minPadj

# geneCount
geneCount


############ $$$$$$$$$$$$$$ ############### 
designAllFactors <- read.delim("./Scripts/TestBasicQC/seqdesign.txt", stringsAsFactors=FALSE)

# Depression VS control within sample NEW

# slice the design and count table and keep less Degraded
moreDegradedDesign <- designAllFactors[designAllFactors$SampleID %in% moreDegraded, ]
moreDegradedNames <- moreDegradedDesign$SampleID
moreDegradedCountTable <- CountTable[,moreDegradedNames]
rownames(moreDegradedCountTable) <- geneIDs

#redo the DEG analysis as before
cds <- newCountDataSet(moreDegradedCountTable, as.factor(moreDegradedDesign$sampleinfo))

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))
cds <- estimateDispersions(cds)
plotDispEsts(cds)

#adjust here#
cds@phenoData@data$condition
res <- nbinomTest(cds, unique(moreDegradedDesign$sampleinfo)[2], unique(moreDegradedDesign$sampleinfo)[1])
# head(res)

completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "moreDegraded_Dep.vs.Cont_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()

minPval <- min(res$pval[!is.na(res$pval)])
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))


## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "moreDegraded_Dep.vs.Cont_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()

completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "moreDegraded_Dep.vs.Cont_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)
moreDegradedRes <- res
# minPval
minPval

# minPadj
minPadj

# geneCount
geneCount

### Draw a compare the foldchange between the less degraded and the more degraded
head(lessDegradedRes)
head(moreDegradedRes)

plot(lessDegradedRes$foldChange, moreDegradedRes$foldChange)
plot(lessDegradedRes$log2FoldChange, moreDegradedRes$log2FoldChange)

##### 
# Select the 1000 most expressed genes, judging by their baseMean
mostExpressedLessDegraded <-  lessDegradedRes[ order(-lessDegradedRes$baseMean), ][1:1000, ]
mostExpressedMoreDegraded <-  moreDegradedRes[ order(-moreDegradedRes$baseMean), ][1:1000, ]
# find genes in common
topExpressedCommonIDs <- mostExpressedMoreDegraded$id[mostExpressedMoreDegraded$id %in% mostExpressedLessDegraded$id]

plot(mostExpressedLessDegraded$foldChange[match( topExpressedCommonIDs, mostExpressedLessDegraded$id)], 
     mostExpressedMoreDegraded$foldChange[match( topExpressedCommonIDs, mostExpressedMoreDegraded$id)])

plot(mostExpressedLessDegraded, mostExpressedMoreDegraded, ylim=c(0,2), xlim=c(0,2))





