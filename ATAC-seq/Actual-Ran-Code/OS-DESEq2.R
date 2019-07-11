#If needed
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#source("https://bioconductor.org/biocLite.R")
#library("genefilter")
#source("https://bioconductor.org/biocLite.R")
#library("pheatmap")
#source("https://bioconductor.org/biocLite.R")
#biocLite("IHW")
#install.packages("whateverismissingtypeitinhereandruntoinstall")
#source("https://bioconductor.org/biocLite.R")
#biocLite("tweeDEseq")

library("DESeq2")
library("kdecopula")
library("foreign")
library("ggplot2")
library("MASS")
library("pheatmap")
library( "gplots" )
library( "RColorBrewer" )

directory <- "/Volumes/Deal-5TB/SUBEXP-ATAC/Rice/UNSCALED-DEseq-CORRECT-NonuniqueAll"
sampleFiles <- grep("counts",list.files(directory),value=TRUE)
sampleFiles
SampleDay <- c("DayOne", "DayTwo", "DayThree", "DayFour", "DayOne", "DayTwo", "DayThree", "DayFour")
sampleCondition <- c("CON", "CON", "CON", "CON", "SUB","SUB", "SUB", "SUB")
Names <- c("OS-CON-1", "OS-CON-2", "OS-CON-3", "OS-CON-4", "OS-SUB-1", "OS-SUB-2", "OS-SUB-3", "OS-SUB-4")
sampleTable <- data.frame(sampleName = Names,
                          fileName = sampleFiles,
                          condition = sampleCondition,
                          Batch = SampleDay)
sampleTable
#VIEWTHETABLE


ddsHTseq <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable,
                                        directory = directory,
                                        design = ~condition) 
colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$condition,
                                      levels = c('CON','SUB'))
dds <- DESeq(ddsHTseq)
res <- results(dds)

resdataL <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdataL)[1] <- 'gene'
head(resdataL)
write.csv(resdataL, file="OS-ATAC-Initial-results.csv")
sum(resdataL$padj < 0.05, na.rm=TRUE )
sum(resdataL$pvalue < 0.05, na.rm=TRUE )

ddsTestTwo <- dds
sizeFactors(ddsTestTwo) <- c(1.227854012254020,	1.218911266300930,	1.001137579065360,	1.117253113693160,	0.676946896941124,	0.891067852160297,	0.790742421354222,	1.076087219941240)
sizeFactors(ddsTestTwo)
ddsTestTwo <- estimateDispersions(ddsTestTwo)
ddsTestTwo <- nbinomWaldTest(ddsTestTwo)
resTestTwo <- results(ddsTestTwo)
plotMA(resTestTwo, ylim = c(-3,3))
dev.copy(png,'OS-MAnorm-Formforallspecies.png')
dev.off()

sum(resTestTwo$padj < 0.05, na.rm=TRUE )
sum(resTestTwo$pvalue < 0.05, na.rm=TRUE )

rld <- rlogTransformation(ddsTestTwo, blind=T)
print(plotPCA(rld, intgroup = c("condition","Batch")))

#STOPPEDHERElasttime


resTestTwoDFr <- merge(as.data.frame(resTestTwo), as.data.frame(counts(ddsTestTwo,normalized=T)), by='row.names',sort=F)
names(resTestTwoDFr)[1] <- 'gene'
head(resTestTwoDFr)
write.csv(resTestTwoDFr, file="OS-resTestTwo-NormalizedResults.csv")



#For All Species - Same Format
with(resTestTwoDFr, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,10), xlim=c(-4,4)))
with(subset(resTestTwoDFr, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resTestTwoDFr, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="gray"))
with(subset(resTestTwoDFr, pvalue<.05 & (log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="purple"))
with(subset(resTestTwoDFr, pvalue<.05 & (log2FoldChange)<(-1)), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))













with(resTestTwoDFr, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", ylim=c(0,20)))
with(subset(resTestTwoDFr, abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(resTestTwoDFr, padj<.1 & abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(resTestTwoDFr, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="gray"))
with(subset(resTestTwoDFr, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))


with(resTestTwoDFr, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,14)))
with(subset(resTestTwoDFr, abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resTestTwoDFr, pvalue<.1 & abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(resTestTwoDFr, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="gray"))
with(subset(resTestTwoDFr, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))







#For All Species - Same Format
with(resTestTwo, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,10), xlim=c(-4,4)))
with(subset(resTestTwo, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resTestTwo, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="gray"))
with(subset(resTestTwo, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.copy(png,'OS-Volcano-Formforallspecies.png')
dev.off()






CHROMS <- read.table("/Users/markobajic/Desktop/RICE-DEseq/CHROMS.txt", header = T, sep = "\t", check.names = F) 
hist(CHROMS$Chromosome, col="blue", breaks=60, freq=NULL, xlim=c(0,13), main = paste("THSs Per Chromosome"), xlab = "Chromosome")

THSsizeAndDistToTSS <- read.table("/Users/markobajic/Desktop/RICE-DEseq/THSsizeAndDistToTSS.txt", header = T, sep = "\t", check.names = F) 
hist(THSsizeAndDistToTSS$THSsize, breaks=600000, freq=NULL, xlim=c(100,500), main = paste("THSs Size"), xlab = "bp")
hist(THSsizeAndDistToTSS$DistanceToTSS, breaks=600000, freq=NULL, xlim=c(-5000,5000), main = paste("Distance To Nearest TSS"), xlab = "bp")


OSUpandDownResults <- read.table("/Users/markobajic/Desktop/RICE-DEseq/OS-UpandDownResults.txt", header = T, sep = "\t", check.names = F) 
hist(OSUpandDownResults$UpInSUBDistToNearestTSS , breaks=4000, freq=NULL, xlim=c(-5000,5000), main = paste("UpInSUB Distance To Nearest TSS"), xlab = "bp")
hist(OSUpandDownResults$DownInSUBDistToNearestTSS , breaks=4000, freq=NULL, xlim=c(-5000,5000), main = paste("DownInSUB Distance To Nearest TSS"), xlab = "bp")
hist(OSUpandDownResults$UpInSUBTHSsize , breaks=600000, freq=NULL, xlim=c(100,500), main = paste("UpInSUB THS Size"), xlab = "bp")
hist(OSUpandDownResults$DownInSUBTHSsize , breaks=600000, freq=NULL, xlim=c(100,500), main = paste("DownInSUB THS Size"), xlab = "bp")





hist(CountsResultsAllFiles$ATCON3, breaks=6000000, freq=NULL, xlim=c(0,200), main = paste("Histogram of ATCON3 counts"), xlab = "counts")
hist(CountsResultsAllFiles$ATSUB3, breaks=6000000, freq=NULL, xlim=c(0,200), main = paste("Histogram of ATSUB3 counts"), xlab = "counts")
hist(CountsResultsAllFiles$ATCON4, breaks=6000000, freq=NULL, xlim=c(0,200), main = paste("Histogram of ATCON4 counts"), xlab = "counts")
hist(CountsResultsAllFiles$ATSUB4, breaks=6000000, freq=NULL, xlim=c(0,200), main = paste("Histogram of ATSUB4 counts"), xlab = "counts")
















#####SUBSET####
directory <- "/Users/markobajic/Desktop/Rice-DEseq/SUBSET"
sampleFiles <- grep("counts",list.files(directory),value=TRUE)
sampleFiles
SampleDay <- c("DayThree", "DayFour", "DayThree", "DayFour")
sampleCondition <- c("CON", "CON", "SUB", "SUB")
Names <- c("SP-SHOOT-CON-3", "SP-SHOOT-CON-4", "SP-SHOOT-SUB-3", "SP-SHOOT-SUB-4")
sampleTable <- data.frame(sampleName = Names,
                          fileName = sampleFiles,
                          condition = sampleCondition,
                          Batch = SampleDay)
sampleTable
#VIEWTHETABLE


ddsHTseq <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable,
                                        directory = directory,
                                        design = ~Batch + condition) 
colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$condition,
                                      levels = c('CON','SUB'))
dds <- DESeq(ddsHTseq)
res <- results(dds)

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=F)), by='row.names', sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
write.csv(resdata, file="OS-3and4-BatchControl.csv")
plotMA(res, ylim = c(-5,5))
sum(res$padj < 0.05, na.rm=TRUE )
sum(res$pvalue < 0.05, na.rm=TRUE )


rld <- varianceStabilizingTransformation(lalla, blind=T)
print(plotPCA(rld, intgroup = c("condition","Batch")))

rld <- rlogTransformation(lalla, blind=T)
print(plotPCA(rld, intgroup = c("condition","Batch")))





ddsTHIS <- dds
sizeFactors(ddsTHIS) <- c(1.211600415,	0.987885118,	1.102463585,	0.879272429,	0.780275046,	1.061842621)
sizeFactors(ddsTHIS)
head(ddsTHIS)

ddsTestOne <- dds
ddsTestOne <- estimateSizeFactors(ddsTestOne)
ddsTestOne <- estimateDispersions(ddsTestOne)
ddsTestOne <- nbinomWaldTest(ddsTestOne)
resTestOne <- results(ddsTestOne)
plotMA(resTestOne, ylim = c(-5,5))
sum(resTestOne$padj < 0.05, na.rm=TRUE )
sum(resTestOne$pvalue < 0.05, na.rm=TRUE )


ddsTestTwo <- dds
sizeFactors(ddsTestTwo) <- c(1.211600415,	0.987885118,	1.102463585,	0.879272429,	0.780275046,	1.061842621)
ddsTestTwo <- estimateDispersions(ddsTestTwo)
ddsTestTwo <- nbinomWaldTest(ddsTestTwo)
resTestTwo <- results(ddsTestTwo)
plotMA(resTestTwo, ylim = c(-5,5))
sum(resTestTwo$padj < 0.05, na.rm=TRUE )
sum(resTestTwo$pvalue < 0.05, na.rm=TRUE )


with(resTestTwo, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", ylim=c(0,6)))
with(subset(resTestTwo, abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(resTestTwo, padj<.1 & abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(resTestTwo, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="gray"))
with(subset(resTestTwo, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))


with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,10)))
with(subset(res, abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, pvalue<.1 & abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="gray"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

















lalla1 <- nbinomWaldTest(dds)
resNoNumbers <- results(lalla1)
lalla2 <- nbinomWaldTest(ddsTHIS)
resWITHNumbers <- results(lalla2)


plotMA(resNoNumbers, ylim = c(-5,5))
sum(resNoNumbers$padj < 0.05, na.rm=TRUE )
sum(resNoNumbers$pvalue < 0.05, na.rm=TRUE )
plotMA(resWITHNumbers, ylim = c(-5,5))
sum(resWITHNumbers$padj < 0.05, na.rm=TRUE )
sum(resWITHNumbers$pvalue < 0.05, na.rm=TRUE )













resdataL <- merge(as.data.frame(resLalla), as.data.frame(counts(lalla,normalized=T)), by='row.names',sort=F)
names(resdataL)[1] <- 'gene'
head(resdataL)
write.csv(resdataL, file="OS-SumAl-DispSizeFDisNegBinomWald.csv")

rld <- varianceStabilizingTransformation(lalla, blind=T)
print(plotPCA(rld, intgroup = c("condition","Batch")))

rld <- rlogTransformation(lalla, blind=T)
print(plotPCA(rld, intgroup = c("condition","Batch")))

















res2 <- results(lalla)
head(lalla)


dds2 <- estimateDispersions(dds2)

lalla <- nbinomWaldTest(dds2)




lalalala <- DESeq(dds2)
res2 <- results(lalalala)
plotMA(res2, ylim = c(-5,5))
sum(res2$padj < 0.05, na.rm=TRUE )
sum(res2$pvalue < 0.05, na.rm=TRUE )
#This did NOT really do anything different, have to figure out just what exactly is going with the adjusting of dds and how to fit my own Size Factors





matrixMB <- read.table("/Users/markobajic/Desktop/RICE-DEseq/NormValues.txt", header = T, sep = "\t", check.names = F) 
dds <- estimateSizeFactors(dds, normMatrix=matrixMB)
head(dds)
head(estimateSizeFactors(dds))


dds <- nbinomLRT(dds)
res2 <- results(dds)

dds <- estimateSizeFactors(dds, normMatrix=m)



sum(res$padj < 0.05, na.rm=TRUE )

sum(res$pvalue < 0.05, na.rm=TRUE )

plotMA(res, ylim = c(-5,5))

dds <- estimateSizeFactors(dds, type = c("poscounts"))
res <- results(dds)
sum(res$padj < 0.05, na.rm=TRUE )
sum(res$pvalue < 0.05, na.rm=TRUE )


dds <- estimateSizeFactors(dds, type = c("poscounts"))
head(sizeFactors(dds))
sizeFactors(dds) <- c(1.211600415,	0.987885118,	1.102463585,	0.879272429,	0.780275046,	1.061842621)
head(sizeFactors(dds))
dds <- sizeFactors(dds)
res <- results(dds)
plotMA(res, ylim = c(-5,5))


res3 <- results(dds)
plotMA(res3, ylim = c(-5,5))
estimateDispersions(dds)
head(dispersions(dds))

nbinomWaldTest(dds, betaPrior = FALSE, modelMatrix = NULL,
               betaTol = 1e-08, maxit = 100, useOptim = TRUE,
               quiet = FALSE, useT = FALSE, df, useQR = TRUE)



library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 500 )

heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))


resultsfile1 <- read.table("/Users/markobajic/Desktop/RICE-DEseq/Rice-Counts-MAR.csv", header = T, sep = ",", check.names = F) 

hist(resultsfile1$ATCON1, breaks=6000000, freq=NULL, xlim=c(0,200), main = paste("Histogram of ATCON1 counts"), xlab = "counts")
hist(resultsfile1$ATSUB2, breaks=6000000, freq=NULL, xlim=c(0,200), main = paste("Histogram of ATSUB2 counts"), xlab = "counts")
hist(resultsfile1$ATCON3, breaks=6000000, freq=NULL, xlim=c(0,200), main = paste("Histogram of ATCON3 counts"), xlab = "counts")
hist(resultsfile1$ATSUB3, breaks=6000000, freq=NULL, xlim=c(0,200), main = paste("Histogram of ATSUB3 counts"), xlab = "counts")
hist(resultsfile1$ATCON4, breaks=6000000, freq=NULL, xlim=c(0,200), main = paste("Histogram of ATCON4 counts"), xlab = "counts")
hist(resultsfile1$ATSUB4, breaks=6000000, freq=NULL, xlim=c(0,200), main = paste("Histogram of ATSUB4 counts"), xlab = "counts")

hist(RandomDevs, breaks=600, freq=NULL, xlim=c(0,80), main = paste("Histogram of RandomDevs"), xlab = "counts")


library("tweeDEseq")

resultsfile1 <- read.table("/Users/markobajic/Desktop/RICE-DEseq/Rice-Counts-MAR.csv", header = T, sep = ",", check.names = F) 
compareCountDist(resultsfile1$"ATCON1", plot=TRUE, xlim=c(10,120))
compareCountDist(RandomDevs, plot=TRUE, xlim=c(10,120))



resultsfile2 <- read.table("/Users/markobajic/Desktop/RICE-DEseq/Rice-Counts-MAR-TestingCutoff.csv", header = T, sep = ",", check.names = F) 
compareCountDist(resultsfile2$"ATCON1", plot=TRUE, xlim=c(10,120))
compareCountDist(RandomDevs, plot=TRUE, xlim=c(10,120))


dPT(resultsfile1$"ATCON1", 44, 2091, 1)
RandomDevs <- rPT(102959, 44, 2091, 1, max = 10*sqrt(44*2091))

goftestfile <- read.table("/Users/markobajic/Desktop/RICE-DEseq/Rice-Counts-MAR.csv", header = T, sep = "\t", check.names = F)
chi2gof <-gofTest(goftestfile, a = 1, mc.cores = 4)
qq <- sort(pchisq(chi2gof, df=1, lower.tail=FALSE))
qqchisq(qq, xlim=c(-1000,1000), ylim=c(-1000,1000))


qqnorm(resultsfile1$ATCON1, main = "ATCON1",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(resultsfile1$ATCON1, datax = FALSE)
qqplot(RandomDevs, resultsfile1$ATCON1)








resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
write.csv(resdata, file="OS-SUBEXP-Results-NoBatchConsideration.csv")

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=F)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
write.csv(resdata, file="OS-SUBEXP-Results-NoBatchConsideration-NOTNORMALIZED.csv")

rld <- varianceStabilizingTransformation(dds, blind=T)
print(plotPCA(rld, intgroup = c("condition","Batch")))

rld <- rlogTransformation(dds, blind=T)
print(plotPCA(rld, intgroup = c("condition","Batch")))

plotSparsity(dds)
plotDispEsts(dds, CV = TRUE)
plotMA(dds)


with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", ylim=c(0,6)))
with(subset(res, abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.1 & abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="gray"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))

with(res_matrix, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", ylim=c(0,6)))
with(subset(res_matrix, abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res_matrix, padj<.1 & abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(res_matrix, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="gray"))
with(subset(res_matrix, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))


with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,8)))
with(subset(res, abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, pvalue<.1 & abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="gray"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

with(res_matrix, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,8)))
with(subset(res_matrix, abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_matrix, pvalue<.1 & abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res_matrix, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="gray"))
with(subset(res_matrix, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))



topVarGenes <- head(order(genefilter::rowVars(assay(rld)), decreasing = TRUE), 20)

mat <- assay(rld) [ topVarGenes, ]
mat <- mat - rowMeans(mat)

df <- as.data.frame(colData(rld))
df <- as.data.frame(colData(rld)[,c("Batch","condition")])

library("pheatmap")
pheatmap(mat, annotation_col=df, cellwidth=40, show_rownames=F, fontsize=12, cluster_cols=F)






dds2 <- nbinomWaldTest(dds)

res2 <- results(dds2)

sum(res2$padj < 0.01, na.rm=TRUE )

resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2,normalized=T)), by='row.names',sort=F)
names(resdata2)[1] <- 'gene'
head(resdata2)
write.csv(resdata2, file="nbinomWaldTest-OS-SUBEXP-Results-NoBatchConsideration.csv")

resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2,normalized=F)), by='row.names',sort=F)
names(resdata2)[1] <- 'gene'
head(resdata2)
write.csv(resdata2, file="nbinomWaldTest-OS-SUBEXP-Results-NoBatchConsideration-NOTNORMALIZED.csv")


rld <- varianceStabilizingTransformation(dds2, blind=T)
print(plotPCA(rld, intgroup = c("condition","Batch")))

rld <- rlogTransformation(dds2, blind=T)
print(plotPCA(rld, intgroup = c("condition","Batch")))

plotSparsity(dds)
plotDispEsts(dds, CV = TRUE)
plotMA(dds)


with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", ylim=c(0,6)))
with(subset(res, abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.1 & abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="gray"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))


with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,6)))
with(subset(res, abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, pvalue<.1 & abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="gray"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))




topVarGenes <- head(order(genefilter::rowVars(assay(rld)), decreasing = TRUE), 20)

mat <- assay(rld) [ topVarGenes, ]
mat <- mat - rowMeans(mat)

df <- as.data.frame(colData(rld))
df <- as.data.frame(colData(rld)[,c("Batch","condition")])

library("pheatmap")
pheatmap(mat, annotation_col=df, cellwidth=40, show_rownames=F, fontsize=12, cluster_cols=F)





