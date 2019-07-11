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

directory <- "/Volumes/Deal-5TB/SUBEXP-ATAC/Tomato2/CountResults/SP-ROOT"
sampleFiles <- grep("counts",list.files(directory),value=TRUE)
sampleFiles
SampleDay <- c("DayTwo", "DayOne", "DayOneB", "DayOne")
sampleCondition <- c("CON", "SUB", "SUB","CON")
Names <- c("SP-ROOT-CON-2", "SP-ROOT-SUB-1", "SP-ROOT-SUB-1B", "SP-ROOT-CON-1")
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

sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm=TRUE )
sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm=TRUE )
sum(res$pvalue < 0.05 & res$log2FoldChange > 1, na.rm=TRUE )
sum(res$pvalue < 0.05 & res$log2FoldChange < -1, na.rm=TRUE )


plotMA(res, ylim = c(-5,3))
dev.copy(png, "SP-ROOT-MAplot-dds.png")
dev.off()

resdataL <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdataL)[1] <- 'gene'
head(resdataL)
write.csv(resdataL, file="SP-ROOT-ATAC-Initial-results.csv")
sum(resdataL$padj < 0.05 & resdataL$log2FoldChange > 1, na.rm=TRUE, na.rm=TRUE )
sum(resdataL$pvalue < 0.05 & resdataL$log2FoldChange < -1, na.rm=TRUE )
sizeFactors(dds)

ddsNORM <- dds
sizeFactors(ddsNORM) <- c(1.262086038332750, 0.723563808153956, 1.007725391790520, 1.006624761722780)
sizeFactors(ddsNORM)      
ddsNORM <- estimateDispersions(ddsNORM)
ddsNORM <- nbinomWaldTest(ddsNORM)
resNORM <- results(ddsNORM)

plotMA(resNORM, ylim = c(-3,3))
dev.copy(png, "SP-ROOT-MAplot-ddsNorm-FORM.png")
dev.off()



sum(resNORM$padj < 0.05 & resNORM$log2FoldChange > 1, na.rm=TRUE )
sum(resNORM$padj < 0.05 & resNORM$log2FoldChange < -1, na.rm=TRUE )
sum(resNORM$pvalue < 0.05 & resNORM$log2FoldChange > 1, na.rm=TRUE )
sum(resNORM$pvalue < 0.05 & resNORM$log2FoldChange < -1, na.rm=TRUE )


rld <- rlogTransformation(ddsNORM, blind=T)
print(plotPCA(rld, intgroup = c("condition","Batch")))
dev.copy(png, "SP-ROOT-PCA-rld-Norm.png")
dev.off()


resNORMdFrame <- merge(as.data.frame(resNORM), as.data.frame(counts(ddsNORM,normalized=T)), by='row.names',sort=F)
names(resNORMdFrame)[1] <- 'gene'
head(resNORMdFrame)
write.csv(resNORMdFrame, file="SP-ROOT-resTestTwo-NormalizedResults.csv")


#For All Species - Same Format
with(resNORMdFrame, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,10), xlim=c(-4,4)))
with(subset(resNORMdFrame, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resNORMdFrame, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="gray"))
with(subset(resNORMdFrame, pvalue<.05 & (log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="purple"))
with(subset(resNORMdFrame, pvalue<.05 & (log2FoldChange)<(-1)), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))












