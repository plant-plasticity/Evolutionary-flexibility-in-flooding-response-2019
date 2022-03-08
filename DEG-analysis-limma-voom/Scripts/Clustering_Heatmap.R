####### DE Gene Clustering ############
## Import limma voom results
limma <-  read.csv("DEG_AllContrasts_171102-subSL.csv", row.names=1, check.names=FALSE)
padj<-0.01            # sets value for cut off
log<-1                # sets value for cut off
sort.col="GeneID"                    # column that will be used to order genes within clusters
col.clust.sort="logFC.SL-HAIRY-TOT" # column on which to order the clusters based on their median value.
col.use=""                          # DataFile columns to use (all="")
## Selection of genes that are differentially expressed, column names can change according to which dataset is selected 
mydata=limma[(limma[,32]<padj&abs(limma[,31])>log)|(limma[,35]<padj&abs(limma[,34])>1)|(limma[,38]<padj&abs(limma[,37])>1),c(37,31,34)]

## Adition of rownames to the selected datasets
rownames(mydata)=unique(as.character(limma[(limma[,32]<padj&abs(limma[,31])>log)|(limma[,35]<padj&abs(limma[,34])>1)|(limma[,38]<padj&abs(limma[,37])>1),c("Symbol.SL-ROOT-TOT")]))
## Library loading
library(cluster)
library(e1071)

## Distance matrix computation
euc <- dist(mydata, method = "euclidean")
## Partition around the mediodids function
pam <-pam(mydata, no.clust, diss = inherits(mydata, "dist"), metric = "euclidean", medoids = NULL, stand = FALSE, cluster.only = FALSE, do.swap = TRUE, pamonce = FALSE, trace.lev = 0)
## Visualization of parameters
summary (pam)
pam$medoids
pam$id.med
pam$clustering
pam$isolation
pam$clusinfo ##gives info on clusters including the no. of genes in each cluster
pam$data
#### User's defined settings
## Selection of number of clusters
no.clust=12
directory=paste0("/MyDirectory/")
date="070719"
name=paste0("dif_SL_root_",date)

setwd(directory)
## Saves clustering file 
write.table(pam$clustering, file=paste0(name,no.clust,"_euc_clust.txt"), sep="\t")
pam1=read.delim(paste0(name,no.clust,"_euc_clust.txt"), sep="\t")
## Merging of clustering files with data
pam2<-merge.data.frame(pam1,mydata,by="row.names",all.x=TRUE, all.y=FALSE)
## Adition of rownames
row.names(pam2)<-pam2[,1]
## Removal of column 1 with names
df2<-pam2[,2:(ncol(mydata)+2)]
## Saves clustering file with data
write.table(df2, file=paste0(name,no.clust,"_euc_clust_val.txt"), sep="\t")
data<-df2
data[1:3,] # check data

## Ser Outputs here
FileTag=paste0(name,"_k",no.clust)  #SPECIFY PREFIX TAG for OUTPUT file names

HeatmapFileFormat=1                 # 1=pdf or 2=tif 
HeatmapFile= paste(FileTag,"_heatmap",sep="") 
## Tag for getting individual files for each cluster
GO_IN=paste(FileTag,"_GO_IN.txt",sep="")

## Obtains median value for cluster
clust.med=cbind(1:no.clust,
                t(sapply(X=1:no.clust,FUN=
                           function(i) {apply(data[which(data[,"x"]==i),2:(ncol(mydata)+1)],
                                              MARGIN=2,FUN=median)})))
clust.med=clust.med[order(clust.med[,col.clust.sort]),]   # sorts the median table and therefore clusters by the selected sample
clust.ord=rev(clust.med[,1])                              # makes vector with descending median of clusters..
clust.med[1:3,]                                           # checkpoint
clust.ord                                                 # checkpoint
## Saves median information
write.table(clust.med,paste(FileTag,"_clust-med.txt",sep=""),sep="\t")

data.ord=matrix(ncol=ncol(data)+1,nrow=0)                 # makes empty matrix
for(i in 1:length(clust.ord)) {
  cluster=cbind(subset(data,
                       data$x==clust.ord[i]),
                "clust.sort"=i)
  data.ord=rbind(data.ord,cluster)}                             #resorts the data and adds a column with the resorted cluster order
data.ord=data.ord[order(data.ord$clust.sort,rownames(data.ord)),]   #sorts by new clust.sort column then by sort col.
#data=data[order(data$cluster,data[,sort.col]),]                    # sorting: number "2" sets the column to sort by after the cluster # if skipping the cluster sorting use these lines
## Saves dataset
write.table(clust.ord,paste(FileTag,"_reorder.txt",sep=""), sep="\t")
data.ord[1:3,]                                            # checkpoint 

data=as.matrix(data.ord[,-c(1,50)])
data[1:3,]

## load libraries 
library(gplots)
library(RColorBrewer)




##Set parameters here
sort.col="GeneID"                    # column that will be used to order genes within clusters
clustGO=read.delim(file=paste0(name,no.clust,"_euc_clust.txt"), sep="\t",row.names=NULL)
colnames(clustGO)=c("GeneID",paste0("k",no.clust,"cluster")) 
clustGO[1:3,] # checkpoint
clustGO_IN=clustGO
for(i in 1:length(clust.ord)) clustGO_IN[,paste0("k",no.clust,"cluster")]=replace(clustGO_IN[,paste0("k",no.clust,"cluster")],which(clustGO[,paste0("k",no.clust,"cluster")]==clust.ord[i]),i) #replaces cluster numbers with reordered cluster numbers
clustGO_IN[1:5,]
clustGO_IN=clustGO_IN[order(clustGO_IN[,paste0("k",no.clust,"cluster")]),] # sorts the file by cluster
clustGO_IN[1:5,]
clustGO_IN_sum=table(clustGO_IN[,2])                                   # makes vector with no. genes in each cluster
clustGO_IN[,3]=clustGO_IN_sum[clustGO_IN[,2]]                          # adds column to data with the no. genes in cluster as third column.
clustGO_IN[1:5,]
colnames(clustGO_IN)=c("GeneID",paste0("k",no.clust,"cluster"),"no.members")                   # renames the column names 
clustGO_IN[1:5,]                                                        # shows chunk of the table
write.table(clustGO_IN,file=GO_IN,row.name=F,sep="\t",quote=F)         # writes the cluster reorder sequence to a file

# makes vector with row numbers to specify cluster breaks position with horizontal white lines on the heatmap
clust.div=clustGO_IN_sum[1]
for(i in 2:length(clustGO_IN_sum)) clust.div=c(clust.div,clust.div[i-1]+clustGO_IN_sum[i])
names(clust.div)=paste("c",1:length(clustGO_IN_sum),sep="")
clust.div
#makes a file with the AGI ID list for each cluster
for(i in 1:max(clustGO_IN[,paste0("k",no.clust,"cluster")]))
  write.table(clustGO_IN[which(clustGO_IN[,paste0("k",no.clust,"cluster")]==i),1],paste(FileTag,"_Clust",as.character(i),"GeneID",sep=""),sep="\t",row.names=F,col.names=F)

#hist(data) # set scale of heatmap based on data range distribution

#open .pdf OR .tif as image target
if(HeatmapFileFormat==1)  pdf(paste(HeatmapFile,".pdf",sep=""), width=6.5,height=8.5) #adjust aspect ration h & w
if(HeatmapFileFormat==2) tiff(paste(HeatmapFile,".tif",sep=""), width=6.5,height=8.5,unit="in",res=300)
#if(HeatmapFileFormat==1) pdf(paste(HeatmapFile,"_clustMed.pdf",sep=""),paper="letter",width=6.5,height=3.5) else tiff(paste(HeatmapFile,".tif",sep=""),width=6.5,height=8.5,unit="in",res=300)

library(RColorBrewer)
#color palette
my_palette <- colorRampPalette(c("#00441b","#1b7837", "#5aae61", "#a6dba0", "#d9f0d3", "#f7f7f7", "#e7d4e8","#c2a5cf",  "#9970ab", "#762a83", "#40004b"))(n = 299)

# original data without gene IDs
data=data[,1:ncol(mydata)]
colnames(data)

## makes the heatmap
pdf(paste0(name,no.clust,"_pam_heatmap.pdf"))
heatmap.2(data,
          #clust.med[rev(1:nrow(clust.med)),2:ncol(clust.med)], # use this to plot collapsed clusters
          col=my_palette,                                       # this sets the color palette made above
          Colv=F,Rowv=F,labRow=F,                               # no dendograms or row labels
          #labCol=fignames[col.use],                            # override column names for figure
          keysize=1,
          cexCol=1,                                             #changes font size of the column labels
          sepcolor="black",
          trace="none",                                         # removes trace lines from heatmap
          density.info="histogram",                             # "histogram" or "none" for the key
          #colsep=c(3,6,11,15,18,21,24,29),   # column indices of vertical white lines 
          rowsep=clust.div,                                     # row indices of horizontal white lines 
          sepwidth=c(0.1,0.1),                             # sets thickness of white lines
          margins=c(10,10),                                       # adjust space for 1-column labels, 2-row labels 
          #lwid=c(0.1,1),                                       # 2x2 relat. widths - Left-key & Rdendo, Right-Cdendo & heatmap
          #lhei=c(0.1,1),                                       # 2x2 relat. heights - Top-key & Cdendo, Bott.-Rdendo & heatmap
          breaks=seq(-6,6,by=12/length(my_palette))             # heat range and bin size - must match no. of colors
)
dev.off()                                                       # close pdf

########### Creates a data.frame to get all genes from each cluster ############
list_clust=data.frame
i=1
x<-paste(FileTag,"_Clust",as.character(i),"GeneID",sep="")
l<-read.delim(x,header=F)
l$clust=paste0("Cluster_",i)

list_clust=l
for(i in 2:no.clust) {
  x<-paste(FileTag,"_Clust",as.character(i),"GeneID",sep="")  # individual files with genes for each cluster
  l<-read.delim(x,header=F)   # read file
  l$clust=paste0("Cluster_",i)  # adition of cluster number
  list_clust=rbind(list_clust,l) # adition of rows for each set
}
colnames(list_clust)=c("gene","comparison") # Add colnames
list_clust[,3:(ncol(mydata)+2)]=mydata[match(list_clust$gene,rownames(mydata)),] # incorporates data of log FC values

## Read Annotation files
anno_Sly=read.csv("/Users/mauricio/Documents/Mauricio JBS postdoc/PUBLICATIONS/SUB_Paper/SUB-paper-tomato/SolycAnnotationsITAG310.csv",header=T)
## Add description column to the dataset (numbers takes into account GeneID length)
list_clust$Description=anno_Sly$Description[match(substr(list_clust$gene,1,14),substr(anno_Sly$Gene,1,14))]

## Read Gene families file
phyto=read.RDS("/phyto_add3_1.RDS")
## Add Gene Family column to the dataset
list_clust$Cluster_ID=phyto$Cluster.ID[match(substr(list_clust$gene,1,14),substr(phyto$Gene.Name,1,14))]

## Save data to a file
write.table(list_clust,paste0(name,no.clust,"_pam_val.xls"),sep="\t",col.names = NA,quote=FALSE)


########### Optional #############
############## GO Analysis for each cluster ###################
######## Generation of GO database as indicated in systemPipeR
library("biomaRt")
library(systemPipeR)
listMarts() # To choose BioMart database
m <- useMart("ENSEMBL_MART_PLANT"); listDatasets(m) 
m <- useMart(biomart="plants_mart", dataset="osativa_eg_gene", host="plants.ensembl.org")

listAttributes(m) # Choose data types you want to download
ids<-getBM(attributes=c("go_id","ensembl_gene_id","namespace_1003","external_gene_name","description"), mart=m)
go <- ids[ids[,3]!="",]; ids[,3] <- as.character(ids[,3])
## correction in gene IDs to match uppercase
go$ensembl_gene_id=gsub("g","G",go$ensembl_gene_id)
go$ensembl_gene_id=gsub("s","S",go$ensembl_gene_id)
dir.create("./GO")
write.table(go, "./GO/GOannotationsBiomart_mod_0410.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
### command from systemPipeR
catdb <- makeCATdb(myfile="./GO/GOannotationsBiomart_mod_0410.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
save(catdb, file="./GO/Os_catdb_0410.RData") 
#############

## load data base containing GO data
load("/MyDirectory/GO/catdb_sl_mod.RData")
r<-data.frame
i=1
x<-paste0("/MyDirectory/interaction-limma/dif_hairy_RS_TF_pam_k12_Clust",i,"GeneID") 
## 
l<-read.delim(x,header=F)
colnames(l)<-paste0("Cluster",i)
## Uses function BatchResult from systemPipeR for cluster 1
BatchResult <- GOCluster_Report(catdb=catdb_sl, setlist=as.data.frame(l), method="all", id_type="gene", CLSZ=2, cutoff=0.05, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
r<-BatchResult
no.clust=12
## Uses function BatchResult from systemPipeR for the rest of the clusters
for(i in seq(no.clust)) {
  x<-paste0("/MyDirectory/dif_hairy_RS_TF_pam_k12_Clust",i,"GeneID")
  l<-read.delim(x,header=F)
  colnames(l)<-paste0("Cluster",i)
  BatchResult <- GOCluster_Report(catdb=catdb_sl_mod, setlist=as.data.frame(l), method="all", id_type="gene", CLSZ=2, cutoff=0.05, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
  ## Add rows with new information
  r<-rbind(r,BatchResult)
}
## Saves output file
write.table(r,file=paste0("/MyDirectory/dif_hairy_RS_TF_GOk_ClustGeneIDall.xls"),quote=FALSE, sep="\t") 
