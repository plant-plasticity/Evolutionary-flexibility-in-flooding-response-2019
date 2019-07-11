############# Comparison of regulated genes under submergence ################
## load of library
library (systemPipeR)
## Set up cutoff value to call as differentially expressed
padj<-0.05
## Set up log FC value to call as differentially expressed
log<-1
## Loading of limma voom outputs containing logFC and adj pvalues for different species and each comparison
limma_SL <-  read.csv("/MyDirectory/DEG_AllContrasts_171102-subSL.csv", row.names=1, check.names=FALSE)
limma_SP <-  read.csv("/MyDirectory/DEG_AllContrasts_171130-subSP.csv", row.names=1, check.names=FALSE)
limma_MT <-  read.csv("/MyDirectory/DEG_AllContrasts_101917-subMT.csv", row.names=1, check.names=FALSE)
limma_OS <-  read.csv("MyDirectory/DEG_AllContrasts-subOS.csv", row.names=1, check.names=FALSE)
## Extraction of genes observed as differentially expressed in different RNA levels
Degs_HR_any=as.character(limma_SL$`Symbol.SL-HAIRY-TOT`[(limma_SL[,2]<padj&limma_SL[,1]>log)|(limma_SL[,5]<padj&limma_SL[,4]>1)|(limma_SL[,8]<padj&limma_SL[,7]>1)|(limma_SL[,11]<padj&limma_SL[,10]>log)])
Degs_SL_any=as.character(limma_SL$`Symbol.SL-ROOT-TOT`[(limma_SL[,32]<padj&limma_SL[,31]>log)|(limma_SL[,35]<padj&limma_SL[,34]>1)|(limma_SL[,38]<padj&limma_SL[,37]>1)])
Degs_SP_any=as.data.frame(limma_SP$`Symbol.SP-ROOT-TOT`[(limma_SP[,2]<padj&limma_SP[,1]>log)|(limma_SP[,5]<padj&limma_SP[,4]>1)|(limma_SP[,8]<padj&limma_SP[,7]>1)])
## Extraction of genes observed as differentially expressed in different RNA levels
Degs_MT_any=as.character(limma_MT$`Symbol.MT-HAIRY-TOT`[(limma_MT[,2]<padj&limma_MT[,1]>log)|(limma_MT[,5]<padj&limma_MT[,4]>1)|(limma_MT[,8]<padj&limma_MT[,7]>1)|(limma_MT[,11]<padj&limma_MT[,10]>log)])
Degs_Rice_any=as.data.frame(limma_OS$Symbol.TRA_SUB[(limma_OS[,2]<padj&limma_OS[,1]>log)|(limma_OS[,5]<padj&limma_OS[,4]>1)|(limma_OS[,8]<padj&limma_OS[,7]>1)|(limma_OS[,11]<padj&limma_OS[,10]>log)])
## Adition of colnames for OS dataframe
colnames(Degs_Rice_any)="RAP"
## Adition of colnames for SP dataframe
colnames(Degs_SP_any)="Sopen"

setwd("/MyDirectory/Supplemental_venn")

######## Load phytozome cluster list ##########
phyto=readRDS("/phyto_add3_1.RDS")
## Loading of table of conversion between RAP and MSU annotation
annoOS=read.delim("/annotation_Os.xls",header=T)
## Loading of table of orthology between Solanum lycopersicum and Solanum pennellii
Sly_Spe=read.delim("/MyDirectory/310topenn20.txt",header=T)

## Adition of MSU gene ID column to DEGs in rice
Degs_Rice_any$MSU=anno$MSU[match(Degs_Rice_any$RAP,anno$RAP)]
## Adition of Slycopersicum gene ID column to DEGs in Spennellii
Degs_SP_any$Sly=Sly_Spe$ITAG3.10.gene.model[match(substr(Degs_SP_any$Sopen,1,14),substr(Sly_Spe$pennelli.ortholog,1,14))]
## Extraction of phytozome  gene family ID "ClusterID" for each DEG
up_rice_any=phyto$ClusterID[phyto$Genes%in%Degs_Rice_any$MSU]
## Usage of the gene name to identify "ClusterID"
up_rm82_any=phyto$ClusterID[substr(phyto$Genes,1,14)%in%substr(Degs_SL_any,1,14)]
up_hairy_any=phyto$ClusterID[substr(phyto$Genes,1,14)%in%substr(Degs_HR_any,1,14)]
up_rpenn_any=phyto$ClusterID[substr(phyto$Genes,1,14)%in%substr(Degs_SP_any$Sly,1,14)]
up_rmed_any=phyto$ClusterID[phyto$Genes%in%Degs_MT_any]
## Names used for identification in venn diagram
names=c("Os","Mt_HR","Sl_HR","Sl","Sp")
## generation of list containing regulated gene families for each species
clus_up_each_any=list(as.character(up_rice_any),up_rmed_any,up_hairy_any,up_rm82_any,up_rpenn_any)
## addition of names
names(clus_up_each_any)=names

## Date to add to the files
date="070719"
## Start of a pdf with named
pdf(paste0("vennplot_up_regulated_",date,".pdf"))
## Overlapper function from systemPipeR
vennset_up_any <- overLapper(clus_up_each_any,type="vennsets")
## VennPlot function systemPipeR 
vennPlot(list(vennset_up_any), mymain="Roots 2h submergence", mysub=paste0("Clusters_up_regulated"),lcol=c("#b35806","#e08214","#fdb863","#fdb863","#8073ac"),lines=c("#b35806","#e08214","#fdb863","#fdb863","#8073ac"), ccol=c("black"))
## finish plot
dev.off()
## creates an element 
n=vennlist(vennset_up_any)
## Obtain data frame with classified gene families
dl_any <- data.frame(ID = rep(names(n), sapply(n, length)),
                 Obs = unlist(n))
## Addition of gene family annotations
dl_any$Description=phyto$Description[match(dl_any$Obs,phyto$Cluster.ID)]
## Save of results 
saveRDS(dl_any,paste0("Clusters_up_anylevel_",date,"_root_shoot_normalized_together.RDS"))