####### Riboseq datasets ###########
# Data obtained from UCDavis facility HiSeq3000 SR50bp
# 3 Replicates perexperiment

######################## GETTING DATA ###################################
# GNU bash, version 4.2.46(2)-release (x86_64-redhat-linux-gnu)
cd Mydirectory/data/
## Example to get data from server
wget http://slims.bioinformatics.ucdavis.edu/Data/2c0g4epsyi/Unaligned/Project_JSMR_L2_RS_MR/INGOLIA06_S4_L002_R1_001.fastq.gz

##################### TRIM ADAPTER #######################################
## Load R
R
############################
#trimming with fastx_clipper
############################
library(BiocParallel); library(BatchJobs)
ftrim=function(x){
  ##
  setwd("/Mydirectory/data")
  ## load libraries
  library(GenomicFeatures);library(GenomicRanges); library(Rsamtools) # Load required libraries.
  library(BiocParallel);library(rtracklayer);library(systemPipeR);library(ShortRead)
  ## File containing the location of fastq.gz files and a column for names
  targets=read.delim("/Mydirectory/data/targetsRS_Mt0721.txt", header=T,comment.char = "#")
  ## Load fastx_toolkit
  moduleload("fastx_toolkit")
  ## Name for files
  input=targets$FileName
  ## command for file uncompressing fastq.gz files
  gzip=paste0("gunzip ",input[x])
  ## command execution
  system(gzip)
  ## designation of output
  output=paste0(input,".trim")
  ## command for adapter trimming
  fastx_command=paste0("fastx_clipper -a 'CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -v -z -Q33 -i ",gsub(".gz","",input[x])," -o ",output[x])
  ## command execution
  system(fastx_command)
}
## File containing the location of fastq.gz files and a column for names
targets=read.delim("/Mydirectory/data/targetsRS_Mt0721.txt", header=T,comment.char = "#")
## Name for files
input=targets$FileName

## Call to the function to execute in multiple cores
funs <- makeClusterFunctionsSLURM("slurm.tmpl")
resources <- list(walltime="60:00:00", ntasks=length(input), ncpus=length(input), memory="10G") # note this assigns 1Gb of Ram per core. If ncpus is 4 then this will amount to 4Gb total
param <- BatchJobsParam(length(input), resources=resources, cluster.functions=funs)
register(param)
bplapply(seq(along=input), ftrim) # Check status with qstat
## Function ftrim can be called with other loop to trim reads

#### QC #############
cd /Mydirectory/data
library(BiocParallel); library(BatchJobs)
fqc=function(x){
  setwd("/Mydirectory/data")
  library(GenomicFeatures);library(GenomicRanges); library(Rsamtools) # Load required libraries.
  library(BiocParallel);library(rtracklayer);library(systemPipeR);library(ShortRead)
  targets=read.delim("/Mydirectory/data/targetsRS_0920.txt", header=T,comment.char = "#")
  input=targets$FileName
  output=paste0(gsub(".gz","",input),".trim.gz")
  fastx_command=paste0("module load fastqc; fastqc ", output[x])
  system(fastx_command)
}
## Call to the function to execute in multiple cores
funs <- makeClusterFunctionsSLURM("slurm.tmpl")
param <- BatchJobsParam(length(output), resources=list(walltime="20:00:00", nodes="1:ppn=4", memory="10gb"), cluster.functions=funs)
register(param)
bplapply(seq(along=input), fqc)
## Function ftrim can be called with other loop to trim reads

#################################################################
#tophat alignment ####
# #############################
## genome fasta
reference <- "/Mydirectory/data/Oryza_sativa.IRGSP-1.0.30.dna.genome"
## genome annotation
gff3<-"/Mydirectory/data/Oryza_sativa.IRGSP-1.0.30.gff3"
## library loading
library(BiocParallel); library(BatchJobs)
falign=function(x){
  library(GenomicFeatures);library(GenomicRanges); library(Rsamtools) # Load required libraries.
  library(BiocParallel);library(rtracklayer);library(systemPipeR);library(ShortRead);library(GenomicFeatures);library(GenomicRanges); library(Rsamtools)

  ## genome fasta
  reference <- "/Mydirectory/data/Oryza_sativa.IRGSP-1.0.30.dna.genome"
  ## genome annotation
  gff3<-"/Mydirectory/data/Oryza_sativa.IRGSP-1.0.30.gff3"
  ##
  targets=read.delim("/Mydirectory/data/targetsRS_Mt0721.txt", header=T,comment.char = "#")
  ## name of files
  input=targets$FileName
  ## trimmed file location
  output=paste0(gsub(".gz","",input),".trim.gz")
  ## tophat command 
  tophat_command <- paste0("module load bowtie2; module load tophat; tophat -p 4 -g 1 --segment-length 15 -i 30 -I 3000  --library-type fr-secondstrand -G ",gff3," -o ", gsub("/data","/results", output[x]), ".tophat ", reference, " ", output[x])
  ## command execution
  system(tophat_command)
  ## Sort bam
  sortBam(file=paste(gsub("/data","/results", output[x]), ".tophat/accepted_hits.bam", sep=""), destination=paste(gsub("/data","/results", output[x]), ".tophat/accepted_hits", sep="")) 
  ## Index bam
  indexBam(paste(gsub("/data","/results", output[x]), ".tophat/accepted_hits.bam", sep=""))
}
## Call to the function to execute in multiple cores
funs <- makeClusterFunctionsSLURM("/bigdata/serreslab/mreynoso/Big_Exp/slurm.tmpl")
param <- BatchJobsParam(length(output), resources=list(walltime="40:00:00", nodes="1:ppn=4", memory="20gb"), cluster.functions=funs)
register(param)
bplapply(seq(along=input), falign)
## Function falign can be called with other loop to trim reads
