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

directory <- "/Volumes/Deal-5TB/SUBEXP-ATAC/Tomato2/CountResults/SL-ROOT"
sampleFiles <- grep("counts",list.files(directory),value=TRUE)
sampleFiles
SampleDay <- c("DayTwo", "DayTwo", "DayOne", "DayOne")
sampleCondition <- c("CON", "SUB", "CON", "SUB")
Names <- c("SL-ROOT-CON-2", "SL-ROOT-SUB-2", "SL-ROOT-CON-1", "SL-ROOT-SUB-1")
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


plotMA(res, ylim = c(-3,8))
dev.copy(png, "SL-ROOT-MAplot-dds.png")
dev.off()

resdataL <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdataL)[1] <- 'gene'
head(resdataL)
write.csv(resdataL, file="SL-ROOT-ATAC-Initial-results.csv")
sum(resdataL$padj < 0.05 & resdataL$log2FoldChange > 1, na.rm=TRUE, na.rm=TRUE )
sum(resdataL$pvalue < 0.05 & resdataL$log2FoldChange < -1, na.rm=TRUE )
sizeFactors(dds)

ddsNORM <- dds
sizeFactors(ddsNORM) <- c(0.17360258427139, 0.42865181224003, 0.32344876712730, 3.07429683636128)
sizeFactors(ddsNORM)      
ddsNORM <- estimateDispersions(ddsNORM)
ddsNORM <- nbinomWaldTest(ddsNORM)
resNORM <- results(ddsNORM)

plotMA(resNORM, ylim = c(-3,3))
dev.copy(png, "SL-ROOT-MAplot-ddsNorm-FORM.png")
dev.off()



sum(resNORM$padj < 0.05 & resNORM$log2FoldChange > 1, na.rm=TRUE )
sum(resNORM$padj < 0.05 & resNORM$log2FoldChange < -1, na.rm=TRUE )
sum(resNORM$pvalue < 0.05 & resNORM$log2FoldChange > 1, na.rm=TRUE )
sum(resNORM$pvalue < 0.05 & resNORM$log2FoldChange < -1, na.rm=TRUE )


rld <- rlogTransformation(ddsNORM, blind=T)
print(plotPCA(rld, intgroup = c("condition","Batch")))
dev.copy(png, "SL-ROOT-PCA-rld-Norm.png")
dev.off()


resNORMdFrame <- merge(as.data.frame(resNORM), as.data.frame(counts(ddsNORM,normalized=T)), by='row.names',sort=F)
names(resNORMdFrame)[1] <- 'gene'
head(resNORMdFrame)
write.csv(resNORMdFrame, file="SL-ROOT-resTestTwo-NormalizedResults.csv")


#For All Species - Same Format
with(resNORMdFrame, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim=c(0,10), xlim=c(-4,4)))
with(subset(resNORMdFrame, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resNORMdFrame, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="gray"))
with(subset(resNORMdFrame, pvalue<.05 & (log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="purple"))
with(subset(resNORMdFrame, pvalue<.05 & (log2FoldChange)<(-1)), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))











with(resNORMdFrame, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot"))
with(subset(resNORMdFrame, abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(resNORMdFrame, padj<.1 & abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(resNORMdFrame, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="gray"))
with(subset(resNORMdFrame, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))


with(resNORMdFrame, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-4,4), ylim=c(0,10)))
with(subset(resNORMdFrame, abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resNORMdFrame, pvalue<.1 & abs(log2FoldChange)>0.6), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(resNORMdFrame, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="gray"))
with(subset(resNORMdFrame, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.copy(png, "SL-ROOT-VOLCANO-FORM.png")
dev.off()







######UNRELATED
B23 <- read.table("/Volumes/Deal-5TB/SUBEXP-ATAC/Tomato2/CountResults/HAIRY/B23-HAIRY-CON-SE-2-counts.txt", header = F, sep = "\t", check.names = F) 
C22 <- read.table("/Volumes/Deal-5TB/SUBEXP-ATAC/Tomato2/CountResults/HAIRY/C22-HAIRY-SUB-2B-combo-counts.txt", header = F, sep = "\t", check.names = F) 
A22 <- read.table("/Volumes/Deal-5TB/SUBEXP-ATAC/Tomato2/CountResults/HAIRY/A22-HAIRY-SUB-1-combo-counts.txt", header = F, sep = "\t", check.names = F) 
B22 <- read.table("/Volumes/Deal-5TB/SUBEXP-ATAC/Tomato2/CountResults/HAIRY/B22-HAIRY-SUB-2-combo-counts.txt", header = F, sep = "\t", check.names = F) 
C21 <- read.table("/Volumes/Deal-5TB/SUBEXP-ATAC/Tomato2/CountResults/HAIRY/C21-HAIRY-CON-2B-combo-counts.txt", header = F, sep = "\t", check.names = F) 
A21 <- read.table("/Volumes/Deal-5TB/SUBEXP-ATAC/Tomato2/CountResults/HAIRY/A21-HAIRY-CON-1-combo-counts.txt", header = F, sep = "\t", check.names = F) 

qqplot(B23$V2, C21$V2, ylab = "C21 CON counts", xlab = "B23 CON counts")
qqplot(B23$V2, C22$V2, ylab = "C22 SUB counts", xlab = "B23 CON counts")
qqplot(A22$V2, C22$V2, ylab = "C22 SUB counts", xlab = "A22 SUB counts")




