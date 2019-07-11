############ Map to chloroplast and then to mitochondria to filter and mapped then to genome #################
##### Fastq files as obtained from the sequencing facility
## Loading of libraries required throughout the procedure
library(GenomicFeatures);library(GenomicRanges); library(Rsamtools) # Load required libraries.
library(BiocParallel);library(rtracklayer);library(systemPipeR);library(ShortRead);library(GenomicFeatures);library(GenomicRanges); library(Rsamtools)
## Set working directory
setwd("MyDirectory/")
## File containing the location of fastq.gz files and a column for names
targets=read.delim("targetsINRNATOTTRA.txt", header=T,comment.char = "#")
input=targets$FileName
library(BiocParallel); library(BatchJobs)
fmap <- function(x) {
  ## Several packages needed 
  library(GenomicFeatures);library(GenomicRanges); library(Rsamtools) # Load required libraries.
  library(BiocParallel);library(rtracklayer);library(systemPipeR);library(ShortRead);library(GenomicFeatures);library(GenomicRanges); library(Rsamtools)
  library(BiocParallel); library(BatchJobs)
  setwd("MyDirectory/")
  targets=read.delim("targetsINRNATOTTRA.txt", header=T,comment.char = "#")
  input=targets$FileName
  name=targets$SampleName
  ## genome file name
  genome<-"OsChloroplast"
  ## command to obtain for tophat2 and bowtie from fasta file obtained from repository
  build<-paste0("module load bowtie2;
                bowtie2-build ",genome,".fa ",genome)
  ## call to the system where the modules have been installed previously
  system(build)
  ## reference
  reference <- paste0("/MyDirectory/",genome)
  ## Output bam file
  bams=paste0(name,"_",genome,".bam")
  ## tophat mapping command
  tophat_command_ch=paste0("module load samtools; module load bowtie2; module load tophat; tophat -p 8 -g 1 -i 50 -I 3000 -o ",name[x], "_" ,genome," ",reference," ", input[x])
  ## execution in the system
  system(tophat_command_ch)
  ## command to convert unmapped reads to fastq. Usage: bedtools bamtofastq [OPTIONS] -i <BAM> -fq <FASTQ> 
  convert=paste0("module load samtools; module load bedtools; bedtools bamtofastq -i ",name[x], "_" ,genome, "/unmapped.bam -fq ",name[x],"no_Chl.fq")
  ## execution in the system
  system(convert)
  ## fastq files input for mapping to mitochondria
  unalign=paste0(name,"no_Chl.fq")
  ## Mitochondria genome build from fasta file of reference
  genome<-"OsMitochondria"
  ## command to obtain for tophat2 and bowtie from fasta file obtained from repository
  build<-paste0("cd /Mydirectory;module load bowtie2; bowtie2-build ",genome,".fa ",genome)
  ## execution in the system
  system(build)
  reference <- paste0("/Mydirectory/",genome)
  ## output name
  bams=paste0(name,"_",genome,".bam")
  ## tophat mapping command
  tophat_command_mt=paste0("module load samtools; module load bowtie2; module load tophat; tophat -p 8 -g 1 -i 50 -I 3000 -o ",name[x], "_" ,genome," ",reference," ", unalign[x])
  ## execution in the system
  system(tophat_command_mt)
  # command to convert unmapped reads to the mitochondria to fastq.
  convert=paste0("module load samtools; module load bedtools; bedtools bamtofastq -i ",name[x], "_" ,genome, "/unmapped.bam -fq ",name[x],"no_Mt.fq")
  ## execution in the system
  system(convert)
  ##  fastq files input for mapping to the genome
  unalign=paste0(name,"no_Mt.fq")
  ## reference
  reference <- "/Mydirectort/Oryza_sativa.IRGSP-1.0.30.dna.genome"
  ## gff3 genome annotation file
  gff3<-"/Mydirectory/Oryza_sativa.IRGSP-1.0.30.gff3"
  ## command for mapping to the genome
  tophat_command <- paste0("module load bowtie2; module load tophat; tophat -p 8 -g 1 -i 50 -I 3000 -G ",gff3," -o ", name[x], ".tophat2 ", reference, " ", unalign[x])
  ## execution in the system
  system(tophat_command)
  ## sort bam files in R
  sortBam(file=paste(name[x], ".tophat2/accepted_hits.bam", sep=""), destination=paste(name[x], ".tophat2/accepted_hits", sep="")) 
  ## index bam files in R
  indexBam(paste(name[x], ".tophat2/accepted_hits.bam", sep=""))
}
## Call to the function to execute in multiple cores
cluster.functions <- makeClusterFunctionsSLURM("slurm.tmpl")
resources <- list(walltime="20:00:00", ntasks=length(input), ncpus=length(input), memory="20G") # note this assigns 1Gb of Ram per core. If ncpus is 4 then this will amount to 4Gb total
param <- BatchJobsParam(length(input), resources=resources, cluster.functions=cluster.functions)
register(param)
bplapply(seq(along=input), fmap)
## Function fmap can be called with other loop to obtain mapped reads
##################################################



