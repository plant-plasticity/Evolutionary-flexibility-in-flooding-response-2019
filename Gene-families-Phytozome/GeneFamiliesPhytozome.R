######## Gene families from Phytozome ############
## Load tsv file downloaded from phytozome Gene families example Viridiplantae
viridiaplantae=read.delim("/MyDirectory/cluster_members_Viridiplantae_5124.tsv")
## install or update packages
install.packages("dplyr")
install.packages("rlang")
install.packages("tidyr")
install.packages("stringr")
## load libraries
library(rlang)
library(tidyr)
library(dplyr)
library(stringr)
## Dissociate Species and ID into separate rows
viridiaplantae=viridiaplantae %>% 
  mutate(Genes = strsplit(as.character(Genes), " ")) %>% 
  unnest(Genes)
## Separate Species into a different columns using an auxiliary dataframe
split_viridiplantae=as.data.frame(str_split_fixed(viridiaplantae$Genes, "_", 2))
class(split_viridiplantae)
## Add a column to indicate species
viridiaplantae$Species=""
## Replace columns with splitted information
viridiaplantae$Species=split_viridiplantae$V1
viridiaplantae$Genes=split_viridiplantae$V2
#################################################################
####### Addition of genes from newer versions of genomes ########
## Load Orphans from new genome 
## file obtained with blastp to Arabidopsis gene sets with 10-6 Evalue
MW=read.csv("/MyDirectory/Tomato_orphans_arabidopsis_MW.csv",header=F)
## Associate to Arabidopsis cluster via match to gene families
MW$ClusterID=viridiaplantae$ClusterID[match(substr(MW$V2,1,9),viridiaplantae$Genes)]
## Add description to genes
MW$Description=viridiaplantae$Description[match(substr(MW$V2,1,9),viridiaplantae$Genes)]
## Check additions
head(MW[,])
## Add Species
MW$Species="Sly"
## colnames if absent
colnames(MW)[2]="Gene.Name.AT"
colnames(MW)[1]="Genes"
## Selection of genes that were absent in phytozome list
Slyc_add_sel=MW[!substr(MW$Genes,1,14)%in%substr(viridiaplantae$Genes,1,14),c(13,14,1,15)]
## Addition of rows to the list in a new dataframe
phyto12virid=rbind(Slyc_add_sel,viridiaplantae)
head(phyto12virid)
## save as RDS file 
saveRDS(phyto12virid,"phyto12update_add3_1.RDS")
#############################################################
