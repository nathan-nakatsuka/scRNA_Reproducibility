
# Vignette Example:

**Steps:**
<br/>
1) Download the following COVID-19 datasets in Seurat object form from: https://atlas.fredhutch.org/fredhutch/covid/ <br/>
-Wilk, Arunachalam, Lee, Wen<br/>
2) Perform differential expression<br/>
3) Perform SumRank Analyses.<br/>
4) Perform Permutations of case/control status and then differential expression and SumRank analyses of the Permutations.<br/>
5) Use the permutations to calibrate the p-values of the real data.<br/>
6) Plot the results.<br/> 
```
# Read in the following functions that can be found in PseudoBulking.R, SumRank.R, CalibratePValueswithPermutations.R, MakePresenceofDataTable.R, and MakingManhattanPlot.R code: AverageMetaData, PseudobulkSeuratObject_Aggregate, GetCommonGenes, SumRank, MakePresenceofDataTable, PermuteCaseControl, CalibratePValueswithPermutations, MakeManhattanPlot

library(Seurat)
library(DESeq2)
library(SeuratDisk)
library(unifed)

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
# Differential Expression (only doing Monocytes here as an example)
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
SumRank(DatasetNames = COVID_DatasetNames, BroadClusterTypes = BroadClusterTypes_COVID, SuffixofDifferentialExpressionOutput="UpReg", CommonGenes=CommonGenes_COVID, ProportionofTopDatasets=ProportionofDatasetstoUse, PresenceofDataTable=PresenceofDataTable_COVID, directory="/home/mydirectory")

# Do Permutations
# Important Note: The code below is extremely slow due to the long amount of time it takes to run differential expression on the datasets 1000 times.
# If possible, it is ideal if you run this code in parallel (i.e. instead of doing the loop of z=1:1000, run each permutation (or sets of 10-100) independently on a cluster, which will make this much faster).
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

# Do SumRank on the permuted differential expression data.
SumRank(DatasetNames = COVID_DatasetNames, BroadClusterTypes = BroadClusterTypes_COVID,
SuffixofDifferentialExpressionOutput=paste0("Permutation_",as.character(PermutationNumber)),
CommonGenes=CommonGenes_COVID, ProportionofTopDatasets=ProportionofDatasetstoUse, PresenceofDataTable=PresenceofDataTable_COVID, directory=paste0("/home/mydirectory/Permutations"))
for(i in 1:length(Datasets)){
rm(list=paste0("avg_exp_",Datasets[i],"_Permutation",as.character(PermutationNumber)))
}
}

# Combine the permutation results
setwd("/home/mydirectory/Permutations")
TotalNumberofPermutations = 1000
ComparisonTable = data.frame(1:(TotalNumberofPermutations*length(CommonGenes_COVID)))
for(z in 1:TotalNumberofPermutations){
    PVal_DirwinHallTable <- read.table(paste(ClusterofInterest,"_CombinedSignedNegLogPValranksNormalized_DirwinHallPValues_Top_",as.character(ProportionofDatasetstoUse),"_ofDatasets_Permutation_",as.character(z),".txt",sep=""),header=T)
    ComparisonTable[((z-1)*length(CommonGenes_COVID)+1):(z*length(CommonGenes_COVID)),1] = PVal_DirwinHallTable$NegLogPValue
}
setwd("/home/mydirectory")
write.table(ComparisonTable,paste0("/home/mydirectory/PermutationComparisonTable_",ClusterofInterest,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

# Compare permutation results to the results of real data to calibrate the p-values of the real data.
PVal_DirwinHallTable <- read.table(paste(ClusterofInterest,"_CombinedSignedNegLogPValranksNormalized_DirwinHallPVals_Top_",as.character(ProportionofDatasetstoUse),"_ofDatasets_UpReg.txt",sep=""),header=T)
FinalPValues = CalibratePValueswithPermutations(CommonGenes=CommonGenes_COVID, ComparisonTable=ComparisonTable,
PVal_DirwinHallTable=PVal_DirwinHallTable)
write.table(FinalPValues,paste0("/home/mydirectory/FinalPValues_COVID_4Datasets_",ClusterofInterest,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

# Plot Manhattan plot
# Note: this is not as good as the results with 16 datasets but with 4 datasets some of the patterns are still clear.
# Note: the significant genes in this plot were jittered substantially so they would be more clear, but they all had infinite -log10(p-value)s.
pdf(paste0("ManhattanPlot_COVID4Datasets_",ClusterofInterest,".txt"))
MakeManhattanPlot(CalibratedPValuesTable=FinalPValues, OtherNegLogPValueCutoff=3, TopValueCutoff=9,
                  Desiredggtitle="Manhattan Plot of COVID-19 vs. Healthy SumRank Differential Expression in CD4 T Cells",jitter_amount=2.5)
dev.off()
```

![image](https://github.com/user-attachments/assets/4c762f06-5d30-480f-810b-4fc96a687cff)

[ManhattanPlot_COVID4Datasets_CD4_T.pdf](https://github.com/user-attachments/files/17297819/ManhattanPlot_COVID4Datasets_CD4_T.pdf)

