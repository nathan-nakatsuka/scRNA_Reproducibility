
# Vignette Example:

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
# Load in the following functions that can be found in PseudoBulking.R and SumRank.R code: AverageMetaData, PseudobulkSeuratObject_Aggregate, GetCommonGenes, SumRank

library(Seurat)
library(DESeq2)

# Load in datasets
Datasets= c("wilk", "arunachalam", "lee", "wen")
for(i in 1:length(Datasets)){
temp = LoadH5Seurat(paste0(Datasets[i],"_2020_processed.HDF5"))
temp$predicted.celltype.l1_2 <- sub(" ", "_", temp$predicted.celltype.l1)
assign(paste0(Datasets[i],"_Seurat),temp)
}

# Do pseudobulking
for(i in 1:length(Datasets)){
temp = get(paste0(Datasets[i],"_Seurat))
avg_temp = PseudobulkSeuratObject_Aggregate(temp, "predicted.celltype.l1_2")
assign(paste0("avg_exp_",Datasets[i]),avg_temp)
}

# Differential Expression (only doing Monocytes here as an example)
# Note: All of these files will be output to the same directory. User can make another directory and output files there if desired.
ClusterofInterest = "Mono"
for(i in 1:length(Datasets)){
    currentTest <- get(paste("avg_exp_",Datasets[i],sep=""))
    avg_temp <- subset(currentTest, subset=predicted.celltype.l1_2==ClusterofInterest)
    Idents(avg_temp) <- "disease_status_standard"
    temp10 <- FindMarkers(avg_temp, ident.1="COVID-19",ident.2="healthy",test.use="DESeq2",group.by = 'disease_status_standard', logfc.threshold = 0, min.pct=0)
    temp10$Gene = rownames(temp10)
    # Remove all mitochondrial genes.
    temp10 = temp10[!grepl('MT-',temp10$Gene),]
    write.table(temp10, file=paste(Datasets[i],"_",ClusterofInterest,"_DESeq2_Pseudobulk_All.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

# Get PresenceofDataTable


# Get common genes
COVID_DatasetNames = c("wilk","arunachalam","lee","wen")
CommonGenes = GetCommonGenes(COVID_DatasetNames,"Mono",PresenceofDataTable_COVID)


```

