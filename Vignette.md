
# Vignette Example:
<br/>

**Steps:**
<br/>
1) Download the following COVID-19 datasets in Seurat object form from: https://atlas.fredhutch.org/fredhutch/covid/ <br/>
-Wilk, Arunachalam, Lee, Wen<br/>
2) Perform differential expression<br/>
3) Perform SumRank Analyses.<br/>
4) Perform Permutations and SumRank analyses of the Permutations.<br/>
5) Use the permutations to calibrate the p-values of the real data.<br/>
6) Plot the results.<br/> 
```
AverageMetaData <- function(orig, avg, f = NULL) {
  f <- f %||% Idents(object = orig)
  f <- as.character(x = f)
  if (!all(colnames(x = avg) %in% unique(x = f))) {
    stop("Not all averaged groups present in original object", call. = FALSE)
  }
  groups <- colnames(x = avg)
  md <- orig[[]]
  for (col in names(x = md)) {
    m <- split(x = md[[col]], f = f)[groups]
    m <- sapply(X = m, FUN = unique, simplify = FALSE, USE.NAMES = TRUE)
    if (!all(vapply(X = m, FUN = length, FUN.VALUE = integer(length = 1L)) == 1L)) {
      warning(
        "Meta data ",
        sQuote(x = col),
        " is not unique across groups",
        call. = FALSE,
        immediate. = TRUE
      )
      next
    }
    message("Adding meta data ", sQuote(x = col), " to averaged object")
    m <- data.frame(unlist(x = m), row.names = names(x = m))
    names(x = m) <- col
    avg[[col]] <- m
  }
  return(avg)
}

# Note, this requires the individual ID to be labeled with the column name "patient".
# Note, this is an aggregate pseudobulk for DESeq2. If you want to get the mean use the PseudobulkSeuratObject_Mean function
PseudobulkSeuratObject_Aggregate <- function(SeuratObject, CellTypeLevel){
   DefaultAssay(SeuratObject) <- "RNA"
   SeuratObject <- NormalizeData(SeuratObject)
   avg_exp_SeuratObject <- Seurat:::PseudobulkExpression(SeuratObject, slot="counts",pb.method='aggregate', group.by = c(CellTypeLevel, "patient"), return.seurat = T)
   #Give back labels to the object
   Idents(SeuratObject) <- CellTypeLevel
   avg_exp_SeuratObject <- AverageMetaData(SeuratObject, avg_exp_SeuratObject, f=paste(as.character(x=Idents(object=SeuratObject)),SeuratObject$patient,sep='_'))
   return(avg_exp_SeuratObject)
}

PseudobulkSeuratObject_Mean <- function(SeuratObject, CellTypeLevel){
   DefaultAssay(SeuratObject) <- "RNA"
   SeuratObject <- NormalizeData(SeuratObject)
   avg_exp_SeuratObject <- AverageExpression(SeuratObject, group.by = c(CellTypeLevel, "patient"), return.seurat = T)
   #Give back labels to the object
   Idents(SeuratObject) <- CellTypeLevel
   avg_exp_SeuratObject <- AverageMetaData(SeuratObject, avg_exp_SeuratObject, f=paste(as.character(x=Idents(object=SeuratObject)),SeuratObject$patient,sep='_'))
   return(avg_exp_SeuratObject)
}

# Load in datasets
Datasets= c("wilk", "arunachalam", "lee", "wen")
for(i in 1:length(Datasets)){
temp = LoadH5Seurat(paste0(Datasets[i],"_2020_processed.HDF5"))
assign(paste0(Datasets[i],"_Seurat),temp)
}

# Do pseudobulking
for(i in 1:length(Datasets)){
temp = get(paste0(Datasets[i],"_Seurat))
avg_temp = PseudobulkSeuratObject_Aggregate(temp, "predicted.celltype.l1")
assign(paste0("avg_exp_",Datasets[i]),avg_temp)
}

# Differential Expression (only doing Monocytes here as an example)
# Note: All of these files will be output to the same directory. User can make another directory and output files there if desired.
ClusterofInterest = "Mono"
for(i in 1:length(Datasets)){
    currentTest <- get(paste("avg_exp_",Datasets[i],sep=""))
    avg_temp <- subset(currentTest, subset=predicted.celltype.l1==ClusterofInterest)
    Idents(avg_temp) <- "disease_status_standard"
    temp10 <- FindMarkers(avg_temp, ident.1="COVID-19",ident.2="healthy",test.use="DESeq2",group.by = 'disease_status_standard', logfc.threshold = 0, min.pct=0)
    temp10$Gene = rownames(temp10)
    # Remove all mitochondrial genes.
    temp10 = temp10[!grepl('MT-',temp10$Gene),]
    write.table(temp10, file=paste(Datasets[i],"_",ClusterofInterest,"_DESeq2_Pseudobulk_All.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}
```

