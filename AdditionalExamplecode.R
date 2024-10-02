library(Seurat)
library(dplyr)
library(Azimuth)
library(SeuratData)
library(presto)
library(car)

#### Merge code
Datasets = c("Mathys","Grubman","Lau","Morabito","Zhou","Leng_EC","OteroGarcia","YangCortex","Gerrits_OTC","Smith_EC","Sadick","Barker","Sayed","SeaAD","Hoffman","Fujita","MathysCell2023")
BroadClusterTypes = c("Oligodendrocytes","Astrocytes","Oligodendrocyte_Precursor_Cell","Glut_ExcitatoryNeuron","Endothelial_Cell","GABA_InhibitoryNeuron","Microglia")
objs <- list()
for(i in 1:length(Datasets)){
	## avg_exp_Datasets[i] is the pseudobulked Seurat object for the particular dataset (e.g. avg_exp_Total above)
	objs[i] <- get(paste("avg_exp_",Datasets[i],sep=""))
}
test.merge <- merge(objs[[1]],objs[2:length(objs)])

BroadClusterTypes = c("Oligodendrocytes","Astrocytes","Oligodendrocyte_Precursor_Cell","Glut_ExcitatoryNeuron","Endothelial_Cell","GABA_InhibitoryNeuron","Microglia")
for(j in 1:length(BroadClusterTypes)){
  #### subset to cluster type.
  test.merge_celltype <- subset(test.merge, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
  test.merge_celltype[["RNA"]]@counts <- as.matrix(test.merge_celltype[["RNA"]]@counts)+1
  test.merge_celltype_DESeq2 <- DESeqDataSetFromMatrix(countData = test.merge_celltype@assays$RNA@counts, colData = test.merge_celltype@meta.data, design = ~Diagnosis+Dataset)
  test.merge_celltype_DESeq2_2 <- DESeq(test.merge_celltype_DESeq2)
  res_ADvsControl <- results(test.merge_celltype_DESeq2_2, contrast=c("Diagnosis","AD","Control"))
  resultsdataframe = data.frame(p_val=res_ADvsControl$pvalue,row.names=rownames(res_ADvsControl),log2fc=res_ADvsControl$log2FoldChange,Gene=rownames(res_ADvsControl))
  write.table(resultsdataframe, file=paste(BroadClusterTypes[j],"_DESeq2_Pseudobulk_MergedDatasets_DEGs_ADvsControl.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}


## SexInteraction Code
MaleDatasets = c("Mathys","Grubman","Lau","Morabito","Zhou","Leng_EC","OteroGarcia","YangCortex","Gerrits_OTC","Smith_EC","Sadick","Sayed","SeaAD","Hoffman","Fujita","MathysCell2023")
FemaleDatasets = c("Mathys","Grubman","Lau","Morabito","Zhou","OteroGarcia","Gerrits_OTC","Smith_EC","Sadick","Barker","Sayed","SeaAD","Hoffman","Fujita","MathysCell2023")
MaleFemaleDatasets = MaleDatasets[MaleDatasets %in% FemaleDatasets]
Datasets = MaleFemaleDatasets
for(i in 1:length(Datasets)){
	for(j in 1:length(BroadClusterTypes)){
		if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets[i],][BroadClusterTypes[j]][[1]]==1){
			Dataset_celltype <- subset(currentTest, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
			Dataset_celltype[["RNA"]]@counts <- as.matrix(Dataset_celltype[["RNA"]]@counts)+1
			## Dataset_celltype[["RNA"]]@counts <- round(Dataset_celltype[["RNA"]]@counts) ## This might be necessary if the counts were not working
			Dataset_celltype_DESeq2 <- DESeqDataSetFromMatrix(countData = Dataset_celltype@assays$RNA@counts, colData = Dataset_celltype@meta.data, design = ~SEX+Diagnosis+SEX:Diagnosis)
			Dataset_celltype_DESeq2$Diagnosis <- relevel(Dataset_celltype_DESeq2$Diagnosis,"Control")
			Dataset_celltype_DESeq2$SEX <- relevel(Dataset_celltype_DESeq2$SEX,"M")
		  Dataset_celltype_DESeq2_2 <- DESeq(Dataset_celltype_DESeq2)
		  }
	res_ADvsControl <- results(Dataset_celltype_DESeq2_2, name=c("SEXF.DiagnosisAD"))
	resultsdataframe = data.frame(p_val=res_ADvsControl$pvalue,row.names=rownames(res_ADvsControl),log2fc=res_ADvsControl$log2FoldChange,Gene=rownames(res_ADvsControl),stderror=res_ADvsControl$lfcSE)
	write.table(resultsdataframe, file=paste(Datasets[i],"_",BroadClusterTypes[j],"_AD_DESeq2_SexInteractionwithAD.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}}

## After this follow the vignette examples to get calibrated p-values for the Sex Interaction values (the results above give female specific scores, while male specific scores are simply reversing the magnitudes of the effects above.)


## Sex Composite Score
Datasets = c("Mathys","Grubman","Lau","Morabito","Zhou","Leng_EC","Leng_SFG","OteroGarcia","Yang1","Yang2","Gerrits_OC","Gerrits_OTC","Smith_SSC","Smith_EC","Sadick","Gabitto_SeaAD","Barker","Sayed","Hoffman","Fujita")
for(i in 1:length(Datasets)){
	## avg_exp_Datasets[i] is the pseudobulked Seurat object for the particular dataset (e.g. avg_exp_Total above)
	currentTest <- get(paste("avg_exp_",Datasets[i],sep=""))
	### MaleAnalyses
	currentTest_Male <- subset(currentTest, subset=SEX=="M")
	for(j in 1:length(BroadClusterTypes)){
		if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets[i],][BroadClusterTypes[j]][[1]]==1){
			avg_temp <- subset(currentTest_Male, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
			### Do Differential expression analysis
			Idents(avg_temp) <- "Diagnosis"
			temp10 <- FindMarkers(avg_temp, ident.1="AD",ident.2="Control",test.use="DESeq2",group.by = 'Diagnosis', logfc.threshold = 0, min.pct=0, min.cells.group=1,pseudocount.use=1)
			temp10$Gene = rownames(temp10)
			temp10 = temp10[!grepl('MT-',temp10$Gene),]
		}
		write.table(temp10, file=paste(Datasets[i],"_",BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_ADvsControlinM.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	}
	### FemaleAnalyses
	currentTest_Female <- subset(currentTest, subset=SEX=="F")
	for(j in 1:length(BroadClusterTypes)){
		if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets[i],][BroadClusterTypes[j]][[1]]==1){
			avg_temp <- subset(currentTest_Female, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
			### Do Differential expression analysis
			Idents(avg_temp) <- "Diagnosis"
			temp10 <- FindMarkers(avg_temp, ident.1="AD",ident.2="Control",test.use="DESeq2",group.by = 'Diagnosis', logfc.threshold = 0, min.pct=0, min.cells.group=1,pseudocount.use=1)
			temp10$Gene = rownames(temp10)
			temp10 = temp10[!grepl('MT-',temp10$Gene),]
		}
        write.table(temp10, file=paste(Datasets[i],"_",BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_ADvsControlinF.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	}
	### FvM in cases 
	currentTest_AD <- subset(currentTest, subset=Diagnosis=="AD")
	for(j in 1:length(BroadClusterTypes)){
		if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets[i],][BroadClusterTypes[j]][[1]]==1){
			avg_temp <- subset(currentTest_AD, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
            ### Do Differential expression analysis
            Idents(avg_temp) <- "Diagnosis"
            temp10 <- FindMarkers(avg_temp, ident.1="F",ident.2="M",test.use="DESeq2",group.by = 'SEX', logfc.threshold = 0, min.pct=0, min.cells.group=1,pseudocount.use=1)
            temp10$Gene = rownames(temp10)
            temp10 = temp10[!grepl('MT-',temp10$Gene),]
		}
		write.table(temp10, file=paste(Datasets[i],"_",BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_FvsMinCases.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	}
	### FvM in controls
	currentTest_Control <- subset(currentTest, subset=Diagnosis=="Control")
	for(j in 1:length(BroadClusterTypes)){
		if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets[i],][BroadClusterTypes[j]][[1]]==1){
			avg_temp <- subset(currentTest_Control, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
            ### Do Differential expression analysis
            Idents(avg_temp) <- "Diagnosis"
            temp10 <- FindMarkers(avg_temp, ident.1="F",ident.2="M",test.use="DESeq2",group.by = 'SEX', logfc.threshold = 0, min.pct=0, min.cells.group=1,pseudocount.use=1)
            temp10$Gene = rownames(temp10)
            temp10 = temp10[!grepl('MT-',temp10$Gene),]
		}
		write.table(temp10, file=paste(Datasets[i],"_",BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_FvsMinControls.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	}
	### MvsF in cases
	currentTest_AD <- subset(currentTest, subset=Diagnosis=="AD")
	for(j in 1:length(BroadClusterTypes)){
		if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets[i],][BroadClusterTypes[j]][[1]]==1){
            avg_temp <- subset(currentTest_AD, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
            ### Do Differential expression analysis
            Idents(avg_temp) <- "Diagnosis"
            temp10 <- FindMarkers(avg_temp, ident.1="M",ident.2="F",test.use="DESeq2",group.by = 'SEX', logfc.threshold = 0, min.pct=0, min.cells.group=1,pseudocount.use=1)
            temp10$Gene = rownames(temp10)
            temp10 = temp10[!grepl('MT-',temp10$Gene),]
		}
		write.table(temp10, file=paste(Datasets[i],"_",BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_MvsFinCases.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	}
	### MvsF in controls
	currentTest_Control <- subset(currentTest, subset=Diagnosis=="Control")
	for(j in 1:length(BroadClusterTypes)){
		if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets[i],][BroadClusterTypes[j]][[1]]==1){
            avg_temp <- subset(currentTest_Control, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
            ### Do Differential expression analysis
            Idents(avg_temp) <- "Diagnosis"
            temp10 <- FindMarkers(avg_temp, ident.1="M",ident.2="F",test.use="DESeq2",group.by = 'SEX', logfc.threshold = 0, min.pct=0, min.cells.group=1,pseudocount.use=1)
            temp10$Gene = rownames(temp10)
            temp10 = temp10[!grepl('MT-',temp10$Gene),]
          }
          write.table(temp10, file=paste(Datasets[i],"_",BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_MvsFinControls.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
    }
}

## After this step, follow the vignettes above to obtain p-values (and then -log10(p-values)).

## Add the relevant values together to get sex specific scores for M and F.
for(j in 1:length(BroadClusterTypes)){
	PVal_DirwinHallTable <- read.table(paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_ADvsControlinM_UpReg_CombinedSignedNegLogPValranksNormalized_DirwinHallPVals_Top_",PercentageofTopDatasets,"_ofDatasets.txt",sep=""),header=T)
	PVal_DirwinHallTable2 <- read.table(paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_ADvsControlinF_UpReg_CombinedSignedNegLogPValranksNormalized_DirwinHallPVals_Top_",PercentageofTopDatasets,"_ofDatasets.txt",sep=""),header=T)
	PVal_DirwinHallTable3 <- read.table(paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_FvsMinCases_UpReg_CombinedSignedNegLogPValranksNormalized_DirwinHallPVals_Top_",PercentageofTopDatasets,"_ofDatasets.txt",sep=""),header=T)
	PVal_DirwinHallTable4 <- read.table(paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_FvsMinControls_UpReg_CombinedSignedNegLogPValranksNormalized_DirwinHallPVals_Top_",PercentageofTopDatasets,"_ofDatasets.txt",sep=""),header=T)
	PVal_DirwinHallTable5 <- read.table(paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_MvsFinCases_UpReg_CombinedSignedNegLogPValranksNormalized_DirwinHallPVals_",PercentageofTopDatasets,"_ofDatasets.txt",sep=""),header=T)
	PVal_DirwinHallTable6 <- read.table(paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_MvsFinControls_UpReg_CombinedSignedNegLogPValranksNormalized_DirwinHallPVals_",PercentageofTopDatasets,"_ofDatasets.txt",sep=""),header=T)
	ComparisonTable = PVal_DirwinHallTable$NegLogPValue+PVal_DirwinHallTable5$NegLogPValue-PVal_DirwinHallTable2$NegLogPValue-PVal_DirwinHallTable6$NegLogPValue
	ComparisonTable2 = PVal_DirwinHallTable2$NegLogPValue+PVal_DirwinHallTable3$NegLogPValue-PVal_DirwinHallTable$NegLogPValue-PVal_DirwinHallTable4$NegLogPValue
    write.table(ComparisonTable, file=paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_SexSpecificityForM_",PercentageofTopDatasets,"_ofDatasets.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
    write.table(ComparisonTable2, file=paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_SexSpecificityForF_",PercentageofTopDatasets,"_ofDatasets.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

## Follow the vignette examples to calibrate p-values by permutations

