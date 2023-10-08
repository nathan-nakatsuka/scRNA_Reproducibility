
library(Seurat)
library(dplyr)
library(Azimuth)
library(SeuratData)
library(presto)
library(car)


##### Step 1: Obtain cell types by mapping to Azimuth reference

### Get Azimuth human motor cortex
setwd("/published_data/human/Azimuth_humanmotorcortex_reference")
CortexRef = readRDS("ref.Rds")

### Table of all Individuals and Metadata
FileNames = read.table("IDs.txt",header=F)
Individuals = unique(FileNames[,1])
IndividualTable = data.frame(1:length(Individuals))
IndividualTable[,1]=Individuals
Metadata = read.table("MetaData.txt",header=T)

### Create Seurat Object for each individual 10X file
for(i in 1:nrow(IndividualTable)){
  temp = Read10X(data.dir= paste("/datadir/",sep=""), extrastuff=paste(IndividualTable[i,1],".",sep=""))
  temp <- CreateSeuratObject(temp)
  temp$Individual <- IndividualTable[i,1]
  assign(paste(IndividualTable[i,1],"_Seurat",sep=""),temp)
}

### Merge all data and then map them to Azimuth reference (integration beforehand is not necessary if mapping to the same reference)
objs <- list()
for(i in 1:nrow(IndividualTable)){
  objs[i] <- get(paste(IndividualTable[i,1],"_Seurat",sep=""))
}
Total_Seurat <- merge(objs[[1]],objs[2:nrow(IndividualTable)])

####QC (example; note: metrics are chosen based on the individual datasets)
Total_Seurat[["percent.mt"]] <- PercentageFeatureSet(Total_Seurat, pattern = "^MT-",assay="RNA")
Total_Seurat <- subset(Total_Seurat, subset = percent.mt <20 & nCount_RNA<65000 & nCount_RNA>200 & nFeature_RNA<9000 & nFeature_RNA>200)

Total_Seurat <- SCTransform(Total_Seurat, vst.flavor = "v2")
Total_Seurat <- RunPCA(Total_Seurat, npcs = 30)
Total_Seurat <- RunUMAP(Total_Seurat, reduction = "pca", dims = 1:30) 

## Find Anchors
Total_Seurat.anchors <- FindTransferAnchors(reference = CortexRef,
                                             query = Total_Seurat, 
                                             dims = 1:30, 
                                             reference.reduction = "refDR", 
                                             features = rownames(CortexRef[['refDR']]@feature.loadings),
                                             k.filter = NA)
#Map to Azimuth cortex reference
Total_Seurat <- MapQuery(
  anchorset = Total_Seurat.anchors,
  query = Total_Seurat,
  reference = CortexRef,
  refdata = list(
    celltype.l2 = "subclass",
    celltype.l1 = "class"
  ), 
  reduction.model = "refUMAP"
)

Total_Seurat$predicted.celltype.l3 <- Total_Seurat$predicted.celltype.l2
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('Astro') = c('Astrocytes')")
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('Endo') = c('Endothelial_Cell')")
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('Oligo') = c('Oligodendrocytes')")
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('L2/3 IT','L5 IT','L5 ET','L5/6 NP','L6 CT','L6 IT','L6 IT Car3','L6b') = c('Glut_ExcitatoryNeuron')")
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('Vip','Sst','Sst Chodl','Sncg','Pvalb','Lamp5') = c('GABA_InhibitoryNeuron')")
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('Micro-PVM') = c('Microglia')")
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('OPC') = c('Oligodendrocyte_Precursor_Cell')")

Metadata = read.table("MetaData.txt",header=T)
row.names(Metadata)=Metadata$Individual
Metadata2 <- Metadata
Metadata2$Individual <- NULL
md <- Metadata[as.character(x = Total_Seurat$Individual), ]
row.names(x = md) <- colnames(Total_Seurat)

Total_Seurat <- AddMetaData(
  object = Total_Seurat,
  metadata = md
)
Total_Seurat$SEX = Total_Seurat$Sex

saveRDS(Total_Seurat,file="Total_Seurat_filteredmapped_annotated.rds")


#### Step 2: Create pseudo-bulked files for DESeq2

### This function gives back the meta-data to the pseudobulked file.
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

### This is for all analyses except DESeq2
DefaultAssay(Total_Seurat) <- "RNA"
Total_Seurat <- NormalizeData(Total_Seurat)
avg_exp_Total <- AverageExpression(Total_Seurat, group.by = c("predicted.celltype.l3", "Individual"), return.seurat = T)
#Give back labels to the object
Idents(Total_Seurat) <- "predicted.celltype.l3"
avg_exp_Total <- AverageMetaData(Total_Seurat, avg_exp_Total, f=paste(as.character(x=Idents(object=Total_Seurat)),Total_Seurat$Individual,sep='_'))
saveRDS(avg_exp_Total, "avg_exp_Total_Mean_PredictedCellTypeL3.rds")

### This is for DESeq2
DefaultAssay(Total_Seurat) <- "RNA"
Total_Seurat <- NormalizeData(Total_Seurat)
avg_exp_Total <- Seurat:::PseudobulkExpression(Total_Seurat, slot="counts", pb.method='aggregate', group.by = c("predicted.celltype.l3", "Individual"), return.seurat = T)
#Give back labels to the object
Idents(Total_Seurat) <- "predicted.celltype.l3"
avg_exp_Total <- AverageMetaData(Total_Seurat, avg_exp_Total, f=paste(as.character(x=Idents(object=Total_Seurat)),Total_Seurat$Individual,sep='_'))
saveRDS(avg_exp_Total, "avg_exp_Total_DESeq2_Aggregate_PredictedCellTypeL3.rds")


#### Step 3. Get differential expression.
Datasets = c("Mathys","Grubman","Lau","Morabito","Zhou","Leng_EC","Leng_SFG","OteroGarcia","Yang1","Yang2","Gerrits_OC","Gerrits_OTC","Smith_SSC","Smith_EC","Sadick","Gabitto_SeaAD","Barker","Sayed","Hoffman","Fujita")
BroadClusterTypes = c("Oligodendrocytes","Astrocytes","Oligodendrocyte_Precursor_Cell","Glut_ExcitatoryNeuron","Endothelial_Cell","GABA_InhibitoryNeuron","Microglia")
for(i in 1:length(Datasets)){
	## avg_exp_Datasets[i] is the pseudobulked Seurat object for the particular dataset (e.g. avg_exp_Total above)
	currentTest <- get(paste("avg_exp_",Datasets[i],sep=""))
	for(j in 1:length(BroadClusterTypes)){
		### PresenceofDataTable is a table that shows whether or not there is data for the particular cell type in that dataset (0=no data; 1=data present)
		if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets[i],][BroadClusterTypes[j]][[1]]==1){
			avg_temp <- subset(currentTest, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
			### Do Differential expression analysis
			Idents(avg_temp) <- "Diagnosis"
			temp10 <- FindMarkers(avg_temp, ident.1="AD",ident.2="Control",test.use="DESeq2",group.by = 'Diagnosis', logfc.threshold = 0, min.pct=0, min.cells.group=1,pseudocount.use=1)
			temp10$Gene = rownames(temp10)
			## Remove all mitochondrial genes
			temp10 = temp10[!grepl('MT-',temp10$Gene),]
    	}
    	write.table(temp10, file=paste(Datasets[i],"_",BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_ADvsControl.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}}



#### Step 4: Rank by Signed NegLog p-value
Datasets = c("Mathys","Grubman","Lau","Morabito","Zhou","Leng_EC","Leng_SFG","OteroGarcia","Yang1","Yang2","Gerrits_OC","Gerrits_OTC","Smith_SSC","Smith_EC","Sadick","Gabitto_SeaAD","Barker","Sayed","Hoffman","Fujita")
BroadClusterTypes = c("Oligodendrocytes","Astrocytes","Oligodendrocyte_Precursor_Cell","Glut_ExcitatoryNeuron","Endothelial_Cell","GABA_InhibitoryNeuron","Microglia")
for(i in 1:length(Datasets)){
	for(j in 1:length(BroadClusterTypes)){
		if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets[i],][BroadClusterTypes[j]][[1]]==1){
			currentTest <- read.table(paste(Datasets[i],"_",BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_ADvsControl.txt",sep=""),header=T)
			## CommonGenes is a vector of genes held in common with all datasets
			currentTest = currentTest[currentTest$Gene %in% CommonGenes,]
			currentTest$NegLogPVal = -log10(currentTest$p_val)
			currentTest$SignedNegLogPVal=currentTest$NegLogPVal
			## Set the SignedNegLogPVal to be negative if the log2FC is negative.
			currentTest$SignedNegLogPVal[currentTest$avg_log2FC<0] <- -currentTest$SignedNegLogPVal[currentTest$avg_log2FC<0]
			currentTest = currentTest[order(currentTest$SignedNegLogPVal,decreasing=TRUE),]
			currentTest$SignedNegLogPVal_rank = rank(-currentTest$SignedNegLogPVal)
			write.table(currentTest, file=paste(Datasets[i],"_",BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_ADvsControl_SignedNegLogPVal_ranked.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
        }}
    }
}}



#### Step 5: Add the signed neglog10p-value ranks of all datasets together
SumRankTable=data.frame(1:length(CommonGenes))
SumRankTable[,1]=CommonGenes[order(as.character(CommonGenes))]
## Loop through each cell type
for(j in 1:length(BroadClusterTypes)){
	for(i in 1:length(Datasets)){
		if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets[i],][BroadClusterTypes[j]][[1]]==1){
			currentTest <- read.table(paste(Datasets[i],"_",BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_ADvsControl_SignedNegLogPVal_ranked.txt",sep=""),header=T)
			## Order all datasets by Gene (so they can be matched properly)
			currentTest = currentTest[order(currentTest$Gene),]
			### Normalize the ranks
			SumRankTable[,i+1]=(currentTest$SignedNegLogPVal_rank-1)/(length(CommonGenes)-1)
		}
		if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets[i],][BroadClusterTypes[j]][[1]]==0){
			SumRankTable[,i+1]=NA
		}
	}
	SumRankTable[,(length(Datasets)+2)]=rowSums(SumRankTable[,2:(length(Datasets)+1)])
	colnames(SumRankTable)=c("Gene",Datasets[1:length(Datasets)],"SumRankAcrossDatasets")
	write.table(SumRankTable, file=paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_AllGenes_CombinedSignedNegLogPValranksNormalized_ADvsControl.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

#### Step 6:  Take the top datasets (choose a percentage) for each gene.
PercentageofTopDatasets = "0.60"
for(j in 1:length(BroadClusterTypes)){
	SumRankTable <- read.table(paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_AllGenes_CombinedSignedNegLogPValranksNormalized_ADvsControl.txt",sep=""),header=T)
	## Get a dataframe without the first and last columns to find the number of relevant Datasets
	test=as.numeric(SumRankTable[1,][,-1])
	test=head(test, -1)
	NumberofRelevantDatasetsForRounding=length(na.exclude(test))
	## Make Smaller tables that fit top percentage of relevant datasets.
	SumRankTable2=SumRankTable
	## Find out how many NA columns there are:
	NumberNAColumns = sum(is.na(test))
	### Sort columns
	sorted_rankings = t(apply(SumRankTable[,2:(ncol(SumRankTable)-1)],MARGIN=1,sort))
	## If there are NA columns, then add them to the end.
	if(NumberNAColumns>0){
		TabletoAdd = data.frame(1:nrow(sorted_rankings))
        TabletoAdd[,1:NumberNAColumns]=NA
        sorted_rankings=cbind(sorted_rankings,TabletoAdd)
	}
	## Get the top percent of datasets by sorting.
	SumRankTable2[,2:(ncol(SumRankTable2)-1)] = sorted_rankings
	Small_SumRankTable=SumRankTable2[,1:(NumberofRelevantDatasetsForRounding*as.numeric(PercentageofTopDatasets)+1)]
	Small_SumRankTable[,((NumberofRelevantDatasetsForRounding*as.numeric(PercentageofTopDatasets))+2)]=rowSums(Small_SumRankTable[,2:((NumberofRelevantDatasetsForRounding*as.numeric(PercentageofTopDatasets))+1)])
	write.table(Small_SumRankTable, file=paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_ADvsControl_CombinedSignedNegLogPValranksNormalized_Top_",PercentageofTopDatasets,"_ofDatasets.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}

#### Step 7: Get p-values for each gene by Irwin Hall distribution
library(unifed)
for(j in 1:length(BroadClusterTypes)){
	SumRankTable <- read.table(paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_ADvsControl_CombinedSignedNegLogPValranksNormalized_Top_",PercentageofTopDatasets,"_ofDatasets.txt",sep=""),header=F)
	PVal_DirwinHallTable=data.frame(1:length(CommonGenes))
	PVal_DirwinHallTable[,1]=CommonGenes[order(as.character(CommonGenes))]
	#### If the sum rank is greater than half, you need to set it to half so it won't be considered significant.
	NumberofRelevantDatasets=ncol(SumRankTable)-2
	SumRankTable[,ncol(SumRankTable)][SumRankTable[,ncol(SumRankTable)] > (NumberofRelevantDatasets/2)] <- NumberofRelevantDatasets/2
	for(i in 1:length(CommonGenes)){
		PVal_DirwinHallTable[i,2]=dirwin.hall(SumRankTable[i,ncol(SumRankTable)],(NumberofRelevantDatasets))
	}
	colnames(PVal_DirwinHallTable)=c("Gene","P_Val")
	PVal_DirwinHallTable$PlotPoint=1:nrow(PVal_DirwinHallTable)
	PVal_DirwinHallTable$NegLogPValue = -log10(PVal_DirwinHallTable$P_Val)
	write.table(PVal_DirwinHallTable, file=paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_ADvsControl_CombinedSignedNegLogPValranksNormalized_DirwinHallPValues_Top_",PercentageofTopDatasets,"_ofDatasets.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

#### Step 8: Calibrate the p-values by Permutation Thresholds
# Final_ComparisonTable = data frame of NegLog10PValues of all Permutations.
for(j in 1:length(BroadClusterTypes)){
Final_ComparisonTable = data.frame(1:1)
PVal_DirwinHallTable <- read.table(paste(BroadClusterTypes[j],"_AD_DESeq2_Pseudobulk_All_ADvsControl_CombinedSignedNegLogPValranksNormalized_DirwinHallPValues_Top_",PercentageofTopDatasets,"_ofDatasets.txt",sep=""),header=T)
FinalTable = data.frame(1:length(CommonGenes))
FinalTable[,1]=CommonGenes
FinalTable[,2:3]=0
ecdf_fun = ecdf(Final_ComparisonTable)
FinalTable[,2]=1-sapply(X = PVal_DirwinHallTable$NegLogPValue, FUN = ecdf_fun)
FinalTable[,3]=p.adjust(FinalTable[,2],method="BH")
colnames(FinalTable)=c("Gene","PVal","PVal_BH")
write.table(FinalTable, file=paste(BroadClusterTypes[j],"AD_DESeq2_Pseudobulk_All_ADvsControl_UpReg_CombinedSignedNegLogPValranksNormalized_DirwinHallPVals_Top_",PercentageofTopDatasets,"_ofDatasets_PermutationThresholdPValues_10000Permutations.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}





#### Merge code
Datasets = c("Mathys","Grubman","Lau","Morabito","Zhou","Leng_EC","Leng_SFG","OteroGarcia","Yang1","Yang2","Gerrits_OC","Gerrits_OTC","Smith_SSC","Smith_EC","Sadick","Gabitto_SeaAD","Barker","Sayed","Hoffman","Fujita")
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
Datasets = c("Mathys","Grubman","Lau","Morabito","Zhou","Leng_EC","Leng_SFG","OteroGarcia","Yang1","Yang2","Gerrits_OC","Gerrits_OTC","Smith_SSC","Smith_EC","Sadick","Gabitto_SeaAD","Barker","Sayed","Hoffman","Fujita")
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

## After this follow steps 4-8 above to get calibrated p-values for the Sex Interaction values (the results above give female specific scores, while male specific scores are simply reversing the magnitudes of the effects above.)


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

## After this step, the results can be plugged into steps 4-7 above to obtain p-values (and then -log10(p-values)).

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

## Follow step 8 above to calibrate p-values by permutations




