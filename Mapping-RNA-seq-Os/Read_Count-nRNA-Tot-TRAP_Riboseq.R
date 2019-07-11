########## Conunt reads over gene features #####################
# Load R
R
## Set working directory
setwd("/Mydirectory/results")
## File containing the location of fastq.gz files and a column for names
targets=read.delim("/Mydirectory/targetsINRNATOTTRA.txt", header=T,comment.char = "#")

fcount <- function(x) {
  ## Loading of libraries
  library(systemPipeR); library(BiocParallel); library(GenomicAlignments);library(GenomicFeatures);library(systemPipeR)
  setwd("/Mydirectory/results")
  targets=read.delim("/Mydirectory/targetsINRNATOTTRA.txt", header=T,comment.char = "#")
  ## Name of the samples
  name=targets$SampleName
  ## location of bam files
  mybams <- paste(name, ".tophat2/accepted_hits.bam", sep="")                                                                                                                                                                    
  ## Name assigned to the bam vector
  names(mybams)=name
  ## Create object with list of bam files
  bfl <- BamFileList(as.character(mybams), yieldSize=50000, index=character())
  ## Assigned to the bam list
  names(bfl)=names(mybams)
  ## Load of Transcript database
  txdb <- loadDb("/Mydirectory/msu.sqlite") 
  ## Functions to extract exons by genes object
  eByg <- exonsBy(txdb, by="gene")
  ## Function to count each bam file over features. Notice data is not strand specific 
  summarizeOverlaps(eByg, bfl[x], mode="Union", ignore.strand=T, inter.feature=FALSE, singleEnd=TRUE)
  ## For Riboseq files instead use the following commented command. This takes into account the strand specificity differences
  # summarizeOverlaps(eByg, bfl[x], mode="Union", ignore.strand=FALSE, inter.feature=FALSE, singleEnd=TRUE)
}
## Call to the function to execute in multiple cores
cluster.functions <- makeClusterFunctionsSLURM("slurm.tmpl")
resources <- list(walltime="20:00:00", ntasks=length(input), ncpus=length(input), memory="1G") # note this assigns 1Gb of Ram per core. If ncpus is 4 then this will amount to 4Gb total
param <- BatchJobsParam(length(input), resources=resources, cluster.functions=cluster.functions)
register(param)
counteByg <- bplapply(seq(along=input), fcount) # Check status with qstat
## Function fcount can be called with other loop to obtain mapped reads

## Assign counts to create a new dataframe 
countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
## Assign rownames corresponding to each gene and colnames corresponding to each sample 
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- name
## Export results
write.table(countDFeByg, "countExonsBygenes_INRNA_TOTAL_TRAP_filteredMtCh.xls", quote=FALSE, sep="\t", col.names = NA)
