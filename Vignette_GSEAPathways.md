
# Vignette Example:

**Steps:**
<br/>
1) Download the following COVID-19 datasets in Seurat object form from: https://atlas.fredhutch.org/fredhutch/covid/ <br/>
-Wilk, Arunachalam, Lee, Wen<br/>
2) Perform differential expression<br/>
3) Perform GSEA SumRank Analyses.
4) Perform Permutations of case/control status and then differential expression and GSEA SumRank analyses of the Permutations.<br/>
5) Use the permutations to calibrate the p-values of the real data.<br/>
```
# Read in the following functions that can be found in PseudoBulking.R, SumRank_GSEAPathways.R, CalibratePValueswithPermutations.R, and MakePresenceofDataTable.R code: AverageMetaData, PseudobulkSeuratObject_Aggregate, GetCommonGenes, SumRank_GSEAPathways, MakePresenceofDataTable, PermuteCaseControl, CalibratePValueswithPermutations

library(Seurat)
library(DESeq2)
library(SeuratDisk)
library(unifed)
library(qvalue)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot) #this package contains visualizations like barplot and dotplot
library(msigdbr)

# Load in datasets
Datasets= c("wilk", "arunachalam", "lee", "wen")
for(i in 1:length(Datasets)){
temp = LoadH5Seurat(paste0(Datasets[i],"_2020_processed.HDF5"))
temp$predicted.celltype.l1_2 <- sub(" ", "_", temp$predicted.celltype.l1)
assign(paste0(Datasets[i],"_Seurat"),temp)
}

# Do pseudobulking
for(i in 1:length(Datasets)){
temp = get(paste0(Datasets[i],"_Seurat"))
avg_temp = PseudobulkSeuratObject_Aggregate(temp, "predicted.celltype.l1_2")
assign(paste0("avg_exp_",Datasets[i]),avg_temp)
}

setwd("/home/mydirectory")
# Differential Expression (only doing CD4_T here as an example)
# Note: This is only for up-regulated genes (for down-regulated genes, just switch the COVID-19 to ident.2 and the healthy to ident.1).
# Note: All of these files will be output to the same directory. User can make another directory and output files there if desired.
ClusterofInterest = "CD4_T"
for(i in 1:length(Datasets)){
    currentTest <- get(paste("avg_exp_",Datasets[i],sep=""))
    avg_temp <- subset(currentTest, subset=predicted.celltype.l1_2==ClusterofInterest)
    Idents(avg_temp) <- "disease_status_standard"
    temp10 <- FindMarkers(avg_temp, ident.1="COVID-19",ident.2="healthy",test.use="DESeq2",group.by = 'disease_status_standard', logfc.threshold = 0, min.pct=0)
    temp10$Gene = rownames(temp10)
    # Remove all mitochondrial genes.
    temp10 = temp10[!grepl('MT-',temp10$Gene),]
    write.table(temp10, file=paste(Datasets[i],"_",ClusterofInterest,"_DifferentialExpression_UpReg.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

# Get PresenceofDataTable
COVID_DatasetNames = c("wilk","arunachalam","lee","wen")
BroadClusterTypes_COVID = unique(avg_exp_wilk$predicted.celltype.l1_2)
PresenceofDataTable_COVID = MakePresenceofDataTable(DatasetNames=COVID_DatasetNames, BroadClusterTypes=BroadClusterTypes_COVID,
CellTypeLevel="predicted.celltype.l1_2")

# Get common genes
CommonGenes_COVID = GetCommonGenes(DatasetNames=COVID_DatasetNames, BroadClusterTypes=BroadClusterTypes_COVID,
                                   PresenceofDataTable=PresenceofDataTable_COVID,CellTypeIndexwithNoMissing = "CD4_T")

ProportionofDatasetstoUse = 1.0
# SumRank
setwd("/home/mydirectory")
CommonGenes=CommonGenes_COVID, ProportionofTopDatasets=ProportionofDatasetstoUse, PresenceofDataTable=PresenceofDataTable_COVID, directory="/home/mydirectory")

# Note: You can change this if you want other biological pathways.
m_df <- msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")
m_df_dataframe <- as.data.frame(m_df)
m_df_dataframe2 = m_df_dataframe[m_df_dataframe$gene_symbol%in% CommonGenes,]
m_df_gene2term = data.frame(TERM=m_df_dataframe2$gs_name,GENE=m_df_dataframe2$gene_symbol)


# This takes ~6.9 minutes.
SumRank_GSEAPathways(DatasetNames = COVID_DatasetNames, BroadClusterTypes = "CD4_T", SuffixofDifferentialExpressionOutput="UpReg", CommonGenes=CommonGenes_COVID, ProportionofTopDatasets=ProportionofDatasetstoUse, PresenceofDataTable=PresenceofDataTable_COVID, directory="/home/mydirectory", m_df_gene2term=m_df_gene2term)

# Do Permutations
# Important Note: The code below is slow (e.g. ~8-12 hours) due to the long amount of time it takes to run differential expression on the datasets 1000 times.
# If possible, it is ideal if you run this code in parallel (i.e. instead of doing the loop of z=1:1000, run each permutation (or sets of 10-100) independently on a cluster, which will make this much faster).
setwd("/home/mydirectory")
FinalPathwayList_Table=read.table("GSEAPathways_FinalPathwayListTable.txt",header=T)

# Note the differential expression does not need to be re-done if it was done on the Vignette.md previously.
setwd("/home/mydirectory/Permutations")
for(z in 1:1000){
PermutationNumber=z
for(i in 1:length(Datasets)){
    currentTest <- get(paste("avg_exp_",Datasets[i],sep=""))
    avg_temp <- subset(currentTest, subset=predicted.celltype.l1_2==ClusterofInterest)
    avg_temp <- PermuteCaseControl(DatasetName=Datasets[i], avg_exp_Dataset=avg_temp, PresenceofDataTable=PresenceofDataTable_COVID,
CellTypeLevel="predicted.celltype.l1_2", BroadClusterTypes="CD4_T", CaseName="COVID-19", ControlName="healthy")
    assign(paste0("avg_exp_",Datasets[i],"_Permutation",as.character(PermutationNumber)),avg_temp)
}

# Do Differential expression on permutations
for(i in 1:length(Datasets)){
    currentTest <- get(paste("avg_exp_",Datasets[i],"_Permutation",as.character(PermutationNumber),sep=""))
    avg_temp <- subset(currentTest, subset=predicted.celltype.l1_2==ClusterofInterest)
    Idents(avg_temp) <- "disease_status_standard_Permuted"
    temp10 <- FindMarkers(avg_temp, ident.1="COVID-19",ident.2="healthy",test.use="DESeq2",group.by = 'disease_status_standard_Permuted', logfc.threshold = 0, min.pct=0)
    temp10$Gene = rownames(temp10)
    # Remove all mitochondrial genes.
    temp10 = temp10[!grepl('MT-',temp10$Gene),]
    write.table(temp10, file=paste(Datasets[i],"_",ClusterofInterest,"_DifferentialExpression_Permutation_",as.character(PermutationNumber),".txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

for(i in 1:length(Datasets)){
rm(list=paste0("avg_exp_",Datasets[i],"_Permutation",as.character(PermutationNumber)))
}

# Do SumRank on the permuted differential expression data.
SumRank_GSEAPathways_Permutations(DatasetNames = COVID_DatasetNames, BroadClusterTypes = "CD4_T",
SuffixofDifferentialExpressionOutput=paste0("Permutation_",as.character(PermutationNumber)),
CommonGenes=CommonGenes_COVID, ProportionofTopDatasets=ProportionofDatasetstoUse, PresenceofDataTable=PresenceofDataTable_COVID, directory=paste0("/home/mydirectory/Permutations"),m_df_gene2term=m_df_gene2term,FinalPathwayList_Table=FinalPathwayList_Table)
}

# Combine the permutation results
## First make sure all the pathway elements are the same.
z=1
setwd("/home/mydirectory/Permutations")
PVal_DirwinHallTable <- read.table(paste(ClusterofInterest,"_CombinedSignedNegLogPValranksNormalized_DirwinHallPValues_Top_",as.character(ProportionofDatasetstoUse),"_ofDatasets_Permutation_",as.character(z),"_GSEAPathwayRanks.txt",sep=""),header=T)
setwd("/home/mydirectory")
PVal_DirwinHallTable2 <- read.table(paste(ClusterofInterest,"_CombinedSignedNegLogPValranksNormalized_DirwinHallPValues_Top_",as.character(ProportionofDatasetstoUse),"_ofDatasets_UpReg_GSEAPathwayRanks.txt",sep=""),header=T)
common_elements <- Reduce(intersect, list(PVal_DirwinHallTable$Gene, PVal_DirwinHallTable2$Gene, FinalPathwayList_Table[,1]))
FinalPathwayList_Table=data.frame(FinalPathwayList_Table[FinalPathwayList_Table[,1] %in% common_elements,])

setwd("/home/mydirectory/Permutations")
TotalNumberofPermutations = 1000
ComparisonTable = data.frame(1:(TotalNumberofPermutations*nrow(FinalPathwayList_Table)))
for(z in 1:TotalNumberofPermutations){
    PVal_DirwinHallTable <- read.table(paste(ClusterofInterest,"_CombinedSignedNegLogPValranksNormalized_DirwinHallPValues_Top_",as.character(ProportionofDatasetstoUse),"_ofDatasets_Permutation_",as.character(z),"_GSEAPathwayRanks.txt",sep=""),header=T)
  PVal_DirwinHallTable=PVal_DirwinHallTable[PVal_DirwinHallTable$Gene %in% common_elements,]
    ComparisonTable[((z-1)*nrow(FinalPathwayList_Table)+1):(z*nrow(FinalPathwayList_Table)),1] = PVal_DirwinHallTable$NegLogPValue
}
setwd("/home/mydirectory")
write.table(ComparisonTable,paste0("/home/mydirectory/PermutationComparisonTable_",ClusterofInterest,"_GSEAPathwayRanks.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

# Compare permutation results to the results of real data to calibrate the p-values of the real data.
PVal_DirwinHallTable <- read.table(paste(ClusterofInterest,"_CombinedSignedNegLogPValranksNormalized_DirwinHallPValues_Top_",as.character(ProportionofDatasetstoUse),"_ofDatasets_UpReg_GSEAPathwayRanks.txt",sep=""),header=T)
PVal_DirwinHallTable=PVal_DirwinHallTable[PVal_DirwinHallTable$Gene %in% common_elements,]
FinalPValues = CalibratePValueswithPermutations(CommonGenes=FinalPathwayList_Table[,1], ComparisonTable=ComparisonTable,
PVal_DirwinHallTable=PVal_DirwinHallTable)
write.table(FinalPValues,paste0("/home/mydirectory/FinalPValues_COVID_4Datasets_",ClusterofInterest,"_GSEAPathwayRanks.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
```
Top pathways from FinalPValues table.

| Pathway | p-value | p-value_BH |
|:-----------|:------------:|------------:|
| data 1     | data 2       | data 3      |
| data 4     | data 5       | data 6      |
| data 7     | data 8       | data 9      |
| data 1     | data 2       | data 3      |
| data 4     | data 5       | data 6      |
| data 7     | data 8       | data 9      |

```
