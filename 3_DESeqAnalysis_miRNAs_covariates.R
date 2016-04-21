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
df <- read.delim("./Results/miRNA_analysis_sel/CountTable.raw.sel.txt",stringsAsFactors=FALSE)

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
  
completePlotPath= paste0("./Results/miRNA_analysis_covariates/images/MAplot.png")
png(completePlotPath, width=635,height=460)
plotMA(res)
graphics.off()
  
minPadj <- min(res$padj[!is.na(res$padj)])
geneCount <- length(which(res$padj[!is.na(res$padj)]<0.05))

## plot Padj Histogram
completePlotPath <- paste0("./Results/miRNA_analysis_covariates/images/histPval.png")
png(completePlotPath, width=635,height=460)
hist(res$pval,breaks=30)
graphics.off()
  
completeTablePath = paste0("./Results/miRNA_analysis_covariates/DESeqTable.txt")
write.table(res, completeTablePath, sep="\t", row.names = FALSE, quote=FALSE)



