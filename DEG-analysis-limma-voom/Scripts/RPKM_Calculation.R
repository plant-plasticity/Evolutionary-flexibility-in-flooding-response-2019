### RPKM calculation over limma voom normalized counts 
## Uses part of the following script:

## RNA seq analysis using limma-voom
# jrodriguezm at ucdavis dot edu
# github: rodriguezmDNA
# brady lab

# modifications from github: reynosoma 

######################## User defined options ########################
######################################################################
# Main directory:
setwd("/Users/mauricio/Documents/Mauricio JBS postdoc/PUBLICATIONS/SUB_Paper/SUB-paper-tomato/SUB-SL-interaction-ALL/") #Full path of the working directory. 
#It must contain: 
# ** A directory named "Counts" with a delimited file for raw counts. Rows are genes and columns samples.
# ** A directory named "meta" with a metadata file with information about the samples. Ideally the number of rows in the metadata is the same as in the raw counts.
# ** A directory named Scripts with this script and the 'functions.R' script.


## Metadata options
metaFile <- "SUB-SL-metadata-ALL-Oct.csv" #Name of metadata file in csv format.
doFilter <- T #F
whichFilter <- c("INRNACONSLROOTE2", "INRNACONSLSHOOTE1","INRNACONSLSHOOTE2","INRNASUBHAIRYROOTI6D","INRNASUBSLROOTE1","INRNASUBSLROOTC2","INRNASUBSLSHOOTE2","TRACONHAIRYROOTT3D","TRACONSLROOTH") #If there are libraries that need to be filtered out (Avoid removing columns manually from the raw counts file)

## Counts file name (with extension) 
countsFile <- "SUB-SL-est-counts-ALL-Oct.txt" #Name of metadata file

shortName <- "171102-subSL" #Short name to append at the end of the filenames. If missing it will append the name of the folder where the scripts where run.

## Filter genes with low expression using CPM counts.
filterByCPM <- T #T
CPMcutoff <- 0.5 #2


## pValue (default = 0.05 ) and absolute logFC (default = 2) to color genes on volcano plots
pValCut=0.05  #0.05  
logCut=2 #2

########################
########################

###

library(edgeR)
library(reshape)
library(gplots)
library(RColorBrewer)
library(calibrate)
library(Glimma) 

## Output
outDir = "SUB-SL-interaction-ALL/"
dir.create(outDir, showWarnings=T)

sink('SUB-SL-interaction-ALL/Outputretestbeforeupload.txt')

geneListsDir = "SUB-SL-interaction-ALL/GeneLists"
dir.create(geneListsDir, showWarnings=T)
#
imgDir = "SUB-SL-interaction-ALL/images/"
dir.create(imgDir, showWarnings=T)
## --

if (is.na(shortName)){
  shortName <- basename(getwd())
}

# Load functions
source("Scripts/functions.R")
######## --- --- --- 


## Start of analysis
####################################################################################
####################################################################################
cat("Reading metadata file \n")

meta <- metaDataProcessing(metaFile,doFilter,whichFilter)
head(meta)

#
cat("Reading counts file:",countsFile,"\n")

GeneCounts <- read.delim(paste0("Counts/",countsFile),row.names = 1)
dim(GeneCounts)

## Check that samples in both counts and metadata are the same.
## Use function filterCounts(counts,meta)
tmp <- filterCounts(GeneCounts,meta)
GeneCounts <- tmp[["counts"]]
meta <- tmp[["meta"]]
rm(tmp)
## --


###### Design matrix
## Convert experimental metadata to factors for the design
experimentFactors <- lapply(apply(meta,2,split,""),unlist)
experimentFactors <- as.data.frame(lapply(experimentFactors,as.factor))

cat ("Create the design with these factors:\n")
print(head(experimentFactors))

###  User modified:
####Simplest design taking into account all possible interactions
Groups <- as.factor(paste0(experimentFactors$Sample,experimentFactors$Treatment,experimentFactors$Genotype,experimentFactors$Tissue))
design <- model.matrix(~0+Groups) 
# Example of an interaction
#design <- model.matrix(~0+experimentFactors$Sample*experimentFactors$Treatment) #Sample*Treatment interaction

## Ensures column names are optimal for the contrast design
fixCols <- paste(c("Groups","experimentFactors","\\$","\\:","\\-",
                   colnames(experimentFactors)),sep="",collapse = "|")

colnames(design) <- gsub(fixCols,"",colnames(design))
head(design)


####################################################################################
cat("Removing genes with 0 counts on all conditions \n")
cat("Initial number of genes:",nrow(GeneCounts),"\n")
rmIDX <- which(rowSums(GeneCounts) == 0)
cat("Removing",length(rmIDX),"genes \n")
GeneCounts <- GeneCounts[-rmIDX,]
cat("Remaining number of genes:",nrow(GeneCounts),"\n")

### Use cpms to uncover lowly expressed genes
dge <- DGEList(counts=GeneCounts,remove.zeros = F)


# Filter genes with low CPMs accross replicates 
cat("Replicates of samples range between:", range(table(Groups)),"\n")

#
if (filterByCPM){
  
  sampleMin <- min(table(Groups))
  cat("Filtering reads with low CPMs ( <",CPMcutoff,") in at least",sampleMin,"replicates \n")
  #
  cpm <- cpm(dge)
  keep.exprs <- rowSums(cpm>CPMcutoff)>=sampleMin
  table(keep.exprs)
  
  
  cat("Removing",table(keep.exprs)[1],"genes \n")
  cat("Remaining number of genes:",table(keep.exprs)[2],"\n")
  
  #
  
  y <- dge[keep.exprs, , keep.lib.size = FALSE]
} else {
  cat("Not doing CPM filtering")
}

### Use voom on the dge object with quantile normalization
v <- voom(y, design, plot = TRUE,normalize.method ="quantile")
## Usage of voom normalized object
r=v
## length of the dataset
dim(r)

## 88 is the length of the dataset 
## revertion of log counts
r$E[,1:88]<-2^r$E[,1:88]
## Generalized: indsamp=length(colnames(r$E))
## r$E[,1:indsamp]<-2^r$E[,1:indsamp]
## new data frame to obtain million reads per sample 
m<-r$targets$lib.size/1000000
## Revertion to calculate normalized counts per gene
r$E=t(t(r$E)*m)
## loading of Slycopersicum transcript database
## Creation from GFF file
#txdb <- makeTxDbFromGFF(file="/bigdata/serreslab/mreynoso/RS_med_tom_1cm/data/ITAG3.10_gene_models.gff", format="gff3",
#                        dataSource="ITAR3.1",
#                        organism="Solanum lycopersicum")
txdb=loadDb("Slycopersicum310.sqlite")
## Extraction of cordinates for exons
eByg <- exonsBy(txdb, by="gene")
## subset for expressed only genes 
eByg2=eByg[names(eByg)%in%substr(rownames(r$E),1,16)]
## Calculation of RPKM using systemPipeR function
rpkmDFeByg <- apply(r$E, 2, function(x) returnRPKM(counts=x, ranges=eByg2))
## Export of RPKM normalized dataset for individual samples
write.table(rpkmDFeByg, "RPKM_SL_after_voom_111717.xls", col.names=NA, quote=FALSE, sep="\t")
## Function to calculate Mean of RPKM 
meanExp <- meanNormalizedExpression(rpkmDFeByg,levels(Groups)) 
## Export of RPKM normalized dataset
write.table(meanExp, "RPKM_SL_Mean_after_voom_111717.xls", col.names=NA, quote=FALSE, sep="\t")
