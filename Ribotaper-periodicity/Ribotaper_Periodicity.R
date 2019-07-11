###### Calculation of periodicity over riboseq reads with Ribotaper ##############
R
## Managing gtf file for annotation file of ribotaper example Medicago
gtf=read.delim("/Mydirectory/Medicago_truncatula.MedtrA17_4.0.39.gtf",comment.char = "#",header=F)
## Edition of gene ID names may be needed to match differences between gff3 and gtf
# medicago_clust$gene2=gsub("Medtr","MTR_",medicago_clust$gene)
## Read annotations
anno_mt=read.delim("/Users/mauricio/Documents/Mauricio JBS postdoc/PUBLICATIONS/SUB_Paper/SUB_Paper_Medicago/interaction-limma/SUB-MT-interaction/Clusters/descr_medicago.txt",header=T)
## Select genes that have not a good annotation "Hypothetical protein" or low confidence "LOW QUALITY" or similar
hypo_mt=(anno_mt[anno_mt$Annotation=="hypothetical protein",])
low_qual=anno_sl[substr(anno_sl$Description,1,10)=="LOW QUALIT",]

## Used when available data from TRAP to select well annotated protein coding genes
RPKM_Mt=read.delim("/RPKM_Mt_Mean.xls",header=T,row.names = 1,sep="\t")
## Calculation of quantiles to exclude lowly expressed and very highly expressed to have a representative set of genes
quantile((RPKM_Mt$TRACONMTHAIRY+RPKM_Mt$TRASUBMTHAIRY)/2,c(0.1,0.95))
## Produce a list of TRAP associated transcripts with confident annotations, filtered between 10 and 1000 RPKM 
write.table(gsub (rownames(RPKM_Mt[((RPKM_Mt$TRACONMTHAIRY+RPKM_Mt$TRASUBMTHAIRY)/2)>10&((RPKM_Mt$TRACONMTHAIRY+RPKM_Mt$TRASUBMTHAIRY)/2)<1000&!rownames(RPKM_Mt)%in%hypo_mt$Gene.ID,]),"Medtr","MTR_"),"genes_sel_mt.txt",quote=F,row.names = F,col.names = F)

## Select confident genes for calculation of cordinates via G
# used outside in unix terminal 
cd '/MyDirectory'
## Selection of gene sets in gtf
grep -f genes_sel.txt Medicago_truncatula.MedtrA17_4.0.39.gtf > selected_Medicago_truncatula.MedtrA17_4.0.39.gtf

######## Creating annotation for Medicago truncatula for Ribotaper
cd /data
/scripts/create_annotations_files.bash 'selected_Medicago_truncatula.MedtrA17_4.0.39.gtf' '/data/JCVI.Medtr.v4.20130313.fasta' false false '/data/anno_ribotaper/' '/bigdata/serreslab/jbazin/ribotaper_1.3/Version_1.3/bedtools_dir/' '/scripts/'
########

##### Running Ribotaper over bam files ######
## load module
module load samtools
## go to directory
cd /results/ribotaper/SUB_MT
## Run Ribotaper Example one merged file
/scripts/create_metaplots.bash /Merged_bams/RSSUBHAIRYROOTMT_merged.bam /data/anno_ribotaper/start_stops_FAR.bed SUBMT /bedtools_dir/ /scripts/