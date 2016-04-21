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

design <- read.delim("./ExperimentInfo/SampleInformation/design.txt",stringsAsFactors=FALSE)
#df <- read.delim("./Results/miRNA_analysis_sel/CountTable.raw.sel.txt",stringsAsFactors=FALSE)
df <- read.delim("./Results/miRNA_analysis_sel/CountTable.raw.sel.txt",stringsAsFactors=FALSE)
dfMatureOnly <- df[grepl(pattern = "mature", df$ID), ]

geneIDs <- df[,c(1)]
CountTable <- data.matrix(df[,-c(1)])

# gene IDs and count table for (pre)mature sequneces only
geneIDsMature <-  dfMatureOnly[,c(1)]
CountTableMature <- data.matrix(dfMatureOnly[,-c(1)])

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


# Depression VS control  within males

############ $$$$$$$$$$$$$$ ############### 
designAllFactors <- read.delim("./Scripts/TestBasicQC/seqdesign.txt", stringsAsFactors=FALSE)

# Depression VS control within sample OLD

# slice the design and count table and keep Old samples only
maleSamplesDesign <- designAllFactors[designAllFactors$Gender=="M" , ]
maleSampleNames <- maleSamplesDesign$SampleID
maleSamplesCountTable <- CountTable[,maleSampleNames]
rownames(maleSamplesCountTable) <- geneIDs

#redo the DEG analysis as before
cds <- newCountDataSet(maleSamplesCountTable, as.factor(maleSamplesDesign$sampleinfo))

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))
cds <- estimateDispersions(cds)
plotDispEsts(cds)

#adjust here#
cds@phenoData@data$condition
res <- nbinomTest(cds, unique(maleSamplesDesign$sampleinfo)[2], unique(maleSamplesDesign$sampleinfo)[1])
# head(res)

completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "maleSamples_Dep.vs.Cont_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()

minPval <- min(res$pval[!is.na(res$pval)])
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))


## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "maleSamples_Dep.vs.Cont_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()

completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "maleSamples_Dep.vs.Cont_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)

# minPval
minPval

# minPadj
minPadj

# geneCount
geneCount


# Depression VS control  within females

############ $$$$$$$$$$$$$$ ############### 
designAllFactors <- read.delim("./Scripts/TestBasicQC/seqdesign.txt", stringsAsFactors=FALSE)

# Depression VS control within sample OLD

# slice the design and count table and keep Old samples only
femaleSamplesDesign <- designAllFactors[designAllFactors$Gender=="F" , ]
femaleSampleNames <- femaleSamplesDesign$SampleID
femaleSamplesCountTable <- CountTable[,femaleSampleNames]
rownames(femaleSamplesCountTable) <- geneIDs

#redo the DEG analysis as before
cds <- newCountDataSet(femaleSamplesCountTable, as.factor(femaleSamplesDesign$sampleinfo))

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))
cds <- estimateDispersions(cds)
plotDispEsts(cds)

#adjust here#
cds@phenoData@data$condition
res <- nbinomTest(cds, unique(femaleSamplesDesign$sampleinfo)[2], unique(femaleSamplesDesign$sampleinfo)[1])
# head(res)

completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "femaleSamples_Dep.vs.Cont_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()

minPval <- min(res$pval[!is.na(res$pval)])
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))


## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "femaleSamples_Dep.vs.Cont_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()

completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "femaleSamples_Dep.vs.Cont_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)

# minPval
minPval

# minPadj
minPadj

# geneCount
geneCount


# Depression VS control  within Young

# slice the design and count table and keep Old samples only
youngSamplesDesign <- designAllFactors[designAllFactors$AgeGroup=="Y", ]
youngSampleNames <- youngSamplesDesign$SampleID
youngSamplesCountTable <- CountTable[,youngSampleNames]
rownames(youngSamplesCountTable) <- geneIDs

#redo the DEG analysis as before
cds <- newCountDataSet(youngSamplesCountTable, as.factor(youngSamplesDesign$sampleinfo))

cds <- estimateSizeFactors(cds)
head(counts(cds,normalized=TRUE))
cds <- estimateDispersions(cds)
plotDispEsts(cds)

#adjust here#
cds@phenoData@data$condition
res <- nbinomTest(cds, unique(youngSamplesDesign$sampleinfo)[2], unique(youngSamplesDesign$sampleinfo)[1])
# head(res)

completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "youngPatSamples_Dep.vs.Cont_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()

minPval <- min(res$pval[!is.na(res$pval)])
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))


## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "youngPatSamples_Dep.vs.Cont_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()

completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "youngPatSamples_Dep.vs.Cont_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)

# minPval
minPval

# minPadj
minPadj

# geneCount
geneCount


# Depression VS control  within Old

# slice the design and count table and keep Old samples only
oldSamplesDesign <- designAllFactors[designAllFactors$AgeGroup=="O", ]
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

completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "oldPatSamples_Dep.vs.Cont_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()

minPval <- min(res$pval[!is.na(res$pval)])
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))


## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "oldPatSamples_Dep.vs.Cont_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()

completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "oldPatSamples_Dep.vs.Cont_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)

# minPval
minPval

# minPadj
minPadj

# geneCount
geneCount


############ $$$$$$$$$$$$$$ ############### 
# group one 
lessDegraded = c("S03", "S04", "S05", "S08", "S09", "S10", "S12", "S18", "S19", "S20")
designAllFactors <- read.delim("./Scripts/TestBasicQC/seqdesign.txt", stringsAsFactors=FALSE)

# Depression VS control within sample NEW

# slice the design and count table and keep less Degraded
lessDegradedDesign <- designAllFactors[designAllFactors$SampleID %in% lessDegraded, ]
lessDegradedNames <- lessDegradedDesign$SampleID
lessDegradedCountTable <- CountTableMature[,lessDegradedNames]
rownames(lessDegradedCountTable) <- geneIDsMature

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

completePlotPath= paste0("./Results/miRNA_analysisIoannis/images/", "lessDegraded_", "MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()

minPval <- min(res$pval[!is.na(res$pval)])
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))


## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysisIoannis/images/", "lessDegraded_", "histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()

completeTablePath = paste0("./Results/miRNA_analysisIoannis/", "lessDegraded_", "DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)
lessDegradedRes <- res

# minPval
minPval

# minPadj
minPadj

# geneCount
geneCount



