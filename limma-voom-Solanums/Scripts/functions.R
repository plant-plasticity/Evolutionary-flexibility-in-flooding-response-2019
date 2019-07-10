## 
#jrm
#last update 2017.04.25

getGeneSymbols <- function(genealiasfile){ 
  geneAlias <- read.delim(paste0("meta/",genealiasfile),header = T, fill = NA,stringsAsFactors = F)
  geneAlias <- split(geneAlias,geneAlias$locus_name)
  ## Reduce symbols, skip some redunant names
  AGI2Symbol <- lapply(geneAlias, function(x){
    idx <- grepl("^(?!A[Tt])",x[,"symbol"],perl = T); #Negates genes starting with At or AT
    if (all(idx == F)){
      paste(x[,"symbol"],collapse = "/")
    } else {
      symbols <- paste(x[idx,"symbol"],collapse = "/") }
  })
  # Return the vector
  
  AGI2Symbol <- unlist(AGI2Symbol)
  return(AGI2Symbol)
}



## This chunk will read metadata files and filter out libraries if needed
metaDataProcessing <- function(metaName,doFilter=F,whichFilter=c("")){
  
  # whichFilter # Vector of libraries to be filtered out. eg., c("TRACONrep1","TOTSUBrep2")
  cat("Reading",metaName,"with filtering options:",doFilter,"\n")
  
  if(doFilter){
    cat ("Filtering libraries:", paste0(whichFilter),"\n")
  }

  # Read and sort 
  metaPath <- paste0("meta/",metaName)
  meta <- read.delim(metaPath,stringsAsFactors = F,as.is = T,sep=",")
  
  meta <- meta[order(meta$Genotype),]
  
  # Filter
  #head(meta)
  if (doFilter){
    rmIDX <- which(meta$Name %in% whichFilter)
    meta <- meta[-rmIDX,]
  }
  
  # Removes whitespaces, converts dashes to underscores.
  meta <- apply(meta,c(1,2),function(x){
    gsub(" |-","_",x)
    })
  
  
  return(as.data.frame(meta,stringsAsFactors = F))
}
####################################################################################
####################################################################################


## Check that samples in both counts and metadata are the same.

filterCounts <- function (counts,meta) { 
  ## Keep columns of counts that have a match with the metadata file 
  # In case some libraries were filtered out of the metadata file
  # This ensures that the counts file won't be manipualted by the user, it will ignore the columns from the analysis but th 
  # original file remains intact.
  cat ("Filter Counts function called \n")

  idxKeep <- intersect(meta$Name, colnames(counts))#
  ## This also keeps the order of the columns the same as the order of samples in the metadata file
  counts <- counts[,idxKeep]
  
  
  if(nrow(meta) != ncol(counts)){
    cat("Check that the metadata file contains the same samples as the counts file \n")
    
    cat("!! -- Attempting to reduce metadata file \n")
    excess <- meta$Name[!meta$Name %in% colnames(counts)] 
    cat("Excess samples: \n", paste0(excess,sep="\n"),"\n")
    meta <- meta[!meta$Name %in% excess,]
  } else (cat("OK - Same samples in counts as in metadata file  \n"))
  
  # Re-check
  if(nrow(meta) == ncol(counts)){cat("OK. \n")}
  
  # Return
  return (list("counts"=counts,"meta"=meta))
  
}

## Visuals
testPalette <- function(mypalette){pie(rep(1, length(mypalette)),col = mypalette)}

#### Different colors available:
Colors13 <- c("#bd9a46",RColorBrewer::brewer.pal(12,"Set3"))
#
ColoresPair <- c("darkgreen", "darkblue")
# 
customColors <- sample(
  c(
    "#bd9a46",
    "#2d66a0",	
    "#ba4343",
    "#84ac38",
    "#267501",	
    "#513267",
    "#0f005b",
    "#3fbfcd",
    "#660099",
    "#1466a2",
    "#f49502",
    "#da0000",
    "#e8a3d2",
    "#554aad","#ffa6ad",
    "#dfce9b","#c99d88",
    "#b85706","#e8b93c",
    "#07ac6a","#193438",
    "#32364f","#b666d2")
)


## Asign colors to each factor of a metadata table.
# Uses the custom color vector by default. 
assignColorsToFactors <- function(factorTable,vectorOfColors=customColors){
  colList <- list()
  for (each in colnames(factorTable)){
    uniQ <- levels(factorTable[,each])
    uniQ
    col <- sample(vectorOfColors) #Reshuffles colors 
    col <- col[1:length(uniQ)]
    names(col) <- uniQ
    colList[[each]] <- col[factorTable[,each]]
  }
  ColorTable <- as.data.frame(do.call("cbind",colList),stringsAsFactors = F)
  return(ColorTable)
}


########
#### Produce MDS plots using labels and a factor from the metadata columns.

#For example produceMDS("Name","Treatment") will use Name column as labels and color the samples by Treatment.
# Another argument is gene.selection, which by default is 'common'. Other accepted value is 'pairwise'
produceMDS <- function(labels,byFactor,gene.selection="common",top=50){
  plotMDS(v, labels=meta[,labels], top=top, 
          col=ColorTable[,byFactor], gene.selection=gene.selection,cex=0.6,
          main=byFactor)}


##### Volcano plots 
# Use a list of DE genes (output from topTable) to draw volcano plots.
# Default coloring: pValCut=0.05 and logCut=2. It can be modified
# makeVolcanoPlots(DEList,pValCut=0.01,logCut=1.5)

makeVolcanoPlots <- function(DEList,pValCut=0.05,logCut=2,plotGenes=T){ 
  
  for (contrast in names(DEList)) {
    print (contrast)#}
    SelectColumns <- c(paste(c("logFC","adj.P.Val","Symbol"),contrast,sep="."))
    # --
    DE <- DEList[[contrast]]
    head(DE)
    
    
    ## Temporary subset and change names
    res <- DE[,colnames(DE) %in% SelectColumns]
    colnames(res) <- c('logFC','adjpVal','Symbol')
    #res$Symbol <- gsub("\\|.*$","",res$Symbol) #If symbol has multiple names, use only the first one.
    #head(res)
    
    # Make a basic volcano plot
    with(res, plot(logFC, -log10(adjpVal), pch=".", 
                   main=paste("Volcano plot of",contrast,sep=" "), 
                   ylim=c(0,max(-log10(res$adjpVal))+5),  col="lightgray", 
                   xlim=c(min(res$logFC)-4,max(res$logFC)+4)))
    
    colors <- c("darkslategray4","seagreen4","tomato4")
    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(res, adjpVal<pValCut), points(logFC, -log10(adjpVal), pch="+", col=colors[1]))
    with(subset(res, abs(logFC)>logCut), points(logFC, -log10(adjpVal), pch="o", col=colors[2]))
    with(subset(res, adjpVal<pValCut & abs(logFC)>logCut), points(logFC, -log10(adjpVal), pch="x", col=colors[3]))
    # Label points with the textxy function from the calibrate plot
    if (plotGenes){
      with(subset(res, adjpVal<.01 & abs(logFC)>3.5), 
      textxy(logFC, -log10(adjpVal), labs=Symbol, cex=.5))
    }
    legend("topright", c("AdjpVal < 0.05",
                         "logFC > 2",
                         "pVal&logFC"),
           fill = colors,border = F,
           cex = 0.8,bty="n")
    abline(v = c(2,-2),col="skyblue",lty = 2)
  }
}


##########
condenseListTables <- function(listDFs) {
  cat ("-- make condenseListTables function called \n")
  # First make an empty table
  uniqRows <- Reduce(union,lapply(listDFs,rownames))
  uniqCols <- unlist(sapply(listDFs, colnames))
  zeroTable <- as.data.frame(matrix(0,
                                    nrow = length(uniqRows),
                                    length(uniqCols)),
                             row.names =uniqRows)
  colnames(zeroTable) = uniqCols
  # Then Fill it
  for (each in names(listDFs)){
    cat ("Filling binary table:", each,"\n")
    zeroTable[rownames(listDFs[[each]]),
              colnames(listDFs[[each]])] <- listDFs[[each]]
  }
  #zeroTable <- apply(zeroTable,c(1,2),as.numeric)
  #zeroTable[is.na(zeroTable)] <- 0
  return(zeroTable)
  cat ("-- binary table done \n")
  cat ("\n")
}


### Get mean of normalized expression values between replicates
meanNormalizedExpression <- function(data=normalizedExpression,uniqGroups=levels(Groups)) { 
  ### Create empty matrix
  meanNorm <- matrix(NA,nrow = nrow(data),ncol = length(uniqGroups),dimnames = list(rownames(data),uniqGroups))
  for( each in uniqGroups){ #Fill the matrix
    cat("Mean of", each,"\n")
    meanNorm[,each] <- apply(data[,grep(each,colnames(data))],1,mean)
  }
  return(meanNorm)
}


