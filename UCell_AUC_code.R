### This script gets the UCell scores for each individual dataset.
library(Seurat)
library(dplyr)
library(UCell)


#### Obtaining meta-data
Datasets = c("Mathys","Grubman","Lau","Morabito","Zhou","Leng_EC","Leng_SFG","OteroGarcia","Yang1","Yang2","Gerrits_OC","Gerrits_OTC","Smith_SSC","Smith_EC","Sadick","Gabitto_SeaAD","Barker","Sayed","Hoffman","Fujita")
BroadClusterTypes = c("Oligodendrocytes","Astrocytes","Oligodendrocyte_Precursor_Cell","Glut_ExcitatoryNeuron","Endothelial_Cell","GABA_InhibitoryNeuron","Microglia")
for(i in 1:length(Datasets)){
  PatientMetaData=read.table(paste(Datasets[i],"_Metadata_byIndividual.txt",sep=""),header=T)
  assign(paste("MetaData_",Datasets[i],sep=""),PatientMetaData)
}

#### Seurat pseudo-bulked files by average expression (not aggregate)
for(i in 1:length(Datasets)){
  currentTest = readRDS(paste0("avg_exp_",as.character(Datasets[i]),"_Mean_PredictedCellTypeL3.rds"))
  assign(paste("avg_exp_Mean_PredictedCellTypeL3_",Datasets[i],sep=""),currentTest)
}

### SignificantGenes_Up.txt are the set of genes being tested for UCell scores with each gene associated with a specific cell type (this can also be a single gene)
## Note: for the single gene analyses, the scores were split by cell type and were not combined across cell types.
for(l in 1:length(Datasets)){
  PredictedCellTypeL3_UpReg <- read.table(paste("SignificantGenes_Up.txt",sep=""),header=T)
  PredictedCellTypeL3_UpRegList = split(PredictedCellTypeL3_UpReg$Gene,PredictedCellTypeL3_UpReg$BroadClusterType)
  currentTest <- get(paste("avg_exp_Mean_PredictedCellTypeL3_",Datasets[l],sep=""))
  ## maxRank is set to 16000 to ensure all genes are being tested even if this decreased the power.
  currentTest <- AddModuleScore_UCell(currentTest, features=PredictedCellTypeL3_UpRegList,maxRank=16000)
  UScoreMatrix = data.frame(1:length(unique(currentTest$Individual)))
  patients = unique(currentTest$Individual)
  UScoreMatrix[,1]= patients
  UScoreMatrix[,2:(length(BroadClusterTypes)+1)]=0
  colnames(UScoreMatrix)=c("Individual",BroadClusterTypes)
  for(j in 1:length(BroadClusterTypes)){
    if(PresenceofDataTable[l,(j+1)]==1){
      avg_temp <- subset(currentTest, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
      for(k in 1:length(patients)){
        avg_temp_Ind = avg_temp@meta.data[avg_temp$Individual==patients[k],]
        if(nrow(avg_temp_Ind)>0){
          CurrentName = paste0(BroadClusterTypes[j],"_UCell")
          UScoreMatrix[k,(j+1)]=avg_temp_Ind[CurrentName]
          ### If gene is not present, need to set UCell score to NA
          if(is.na(PredictedCellTypeL3_UpReg[PredictedCellTypeL3_UpReg$BroadClusterType==BroadClusterTypes[j],]$Gene[1])){UScoreMatrix[k,(j+1)]=NA}
        }
        if(nrow(avg_temp_Ind)==0){
          UScoreMatrix[k,(j+1)]=NA
        }
      }
    }
    if(PresenceofDataTable[l,(j+1)]==0){
      UScoreMatrix[,(j+1)]=NA
    }
  }
  	## Add in patient metadata
	PatientMetaData <- get(paste("MetaData_",Datasets[l],sep=""))
	UScoreMatrix_Expanded = cbind(UScoreMatrix,PatientMetaData)
	DatasetSpecificScores=UScoreMatrix_Expanded[c("Individual","Disease_Severity","Diagnosis",BroadClusterTypes)]
	write.table(DatasetSpecificScores, file=paste0("AD_UCellScores_UpRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testing_",Datasets[l],"_Dataset.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
  
  PredictedCellTypeL3_DownReg <- read.table(paste("SignificantGenes_Down.txt",sep=""),header=T)
  PredictedCellTypeL3_DownRegList = split(PredictedCellTypeL3_DownReg$Gene,PredictedCellTypeL3_DownReg$BroadClusterType)
  currentTest <- get(paste("avg_exp_Mean_PredictedCellTypeL3_",Datasets[l],sep=""))
  ## maxRank is set to 16000 to ensure all genes are being tested even if this decreased the power.
  currentTest <- AddModuleScore_UCell(currentTest, features=PredictedCellTypeL3_DownRegList,maxRank=16000)
  UScoreMatrix = data.frame(1:length(unique(currentTest$Individual)))
  patients = unique(currentTest$Individual)
  UScoreMatrix[,1]= patients
  UScoreMatrix[,2:(length(BroadClusterTypes)+1)]=0
  colnames(UScoreMatrix)=c("Individual",BroadClusterTypes)
  for(j in 1:length(BroadClusterTypes)){
    if(PresenceofDataTable[l,(j+1)]==1){
      avg_temp <- subset(currentTest, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
      for(k in 1:length(patients)){
        avg_temp_Ind = avg_temp@meta.data[avg_temp$Individual==patients[k],]
        if(nrow(avg_temp_Ind)>0){
          CurrentName = paste0(BroadClusterTypes[j],"_UCell")
          UScoreMatrix[k,(j+1)]=avg_temp_Ind[CurrentName]
          ### If gene is not present, need to set UCell score to NA
          if(is.na(PredictedCellTypeL3_DownReg[PredictedCellTypeL3_DownReg$BroadClusterType==BroadClusterTypes[j],]$Gene[1])){UScoreMatrix[k,(j+1)]=NA}
        }
        if(nrow(avg_temp_Ind)==0){
          UScoreMatrix[k,(j+1)]=NA
        }
      }
    }
    if(PresenceofDataTable[l,(j+1)]==0){
      UScoreMatrix[,(j+1)]=NA
    }
  }
  	## Add in patient metadata
    PatientMetaData <- get(paste("MetaData_",Datasets[l],sep=""))
    UScoreMatrix_Expanded = cbind(UScoreMatrix,PatientMetaData)
	DatasetSpecificScores=UScoreMatrix_Expanded[c("Individual","Disease_Severity","Diagnosis",BroadClusterTypes)]
	write.table(DatasetSpecificScores, file=paste0("AD_UCellScores_DownRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testing_",Datasets[l],"_Dataset.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

### Normalize all scores.
for(i in 1:length(Datasets)){
  DatasetSpecificScores=read.table(paste0("AD_UCellScores_UpRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testing_",Datasets[i],"_Dataset.txt"),header=T)
  for(j in 1:length(BroadClusterTypes)){
    ClusterMean = mean(DatasetSpecificScores[BroadClusterTypes[j]][[1]], na.rm=TRUE)
    ClusterSum = sum(DatasetSpecificScores[BroadClusterTypes[j]][[1]], na.rm=TRUE)
    ## Only Normalize if there are some scores there for that cell type
    if(ClusterSum>0){
      ClusterMin = min(DatasetSpecificScores[BroadClusterTypes[j]][[1]],na.rm=TRUE)
      ClusterMax = max(DatasetSpecificScores[BroadClusterTypes[j]][[1]],na.rm=TRUE)
      ClusterRange = ClusterMax-ClusterMin
      DatasetSpecificScores[BroadClusterTypes[j]]=(DatasetSpecificScores[BroadClusterTypes[j]]-ClusterMin)/ClusterRange
    }
  }
  write.table(DatasetSpecificScores, file=paste0("AD_UCellScores_UpRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testing_",Datasets[i],"_Dataset_Normalized.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

### Normalize all scores.
for(i in 1:length(Datasets)){
  DatasetSpecificScores=read.table(paste0("AD_UCellScores_DownRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testing_",Datasets[i],"_Dataset.txt"),header=T)
  for(j in 1:length(BroadClusterTypes)){
    ClusterMean = mean(DatasetSpecificScores[BroadClusterTypes[j]][[1]], na.rm=TRUE)
    ClusterSum = sum(DatasetSpecificScores[BroadClusterTypes[j]][[1]], na.rm=TRUE)
    ## Only Normalize if there are some scores there for that cell type
    if(ClusterSum>0){
      ClusterMin = min(DatasetSpecificScores[BroadClusterTypes[j]][[1]],na.rm=TRUE)
      ClusterMax = max(DatasetSpecificScores[BroadClusterTypes[j]][[1]],na.rm=TRUE)
      ClusterRange = ClusterMax-ClusterMin
      DatasetSpecificScores[BroadClusterTypes[j]]=(DatasetSpecificScores[BroadClusterTypes[j]]-ClusterMin)/ClusterRange
    }
  }
  write.table(DatasetSpecificScores, file=paste0("AD_UCellScores_DownRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testing_",Datasets[i],"_Dataset_Normalized.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

###########################
### If there is no score present, then just set it to the mean of all other values (or to 0 if all of the data is not present)
for(i in 1:length(Datasets)){
  DatasetSpecificScores=read.table(paste0("AD_UCellScores_UpRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testing_",Datasets[i],"_Dataset_Normalized.txt"),header=T)
  for(j in 1:length(BroadClusterTypes)){
    Total_Mean = mean(DatasetSpecificScores[BroadClusterTypes[j]][[1]], na.rm=TRUE)
    for(k in 1:nrow(DatasetSpecificScores)){
      if(is.na(DatasetSpecificScores[BroadClusterTypes[j]][[1]][k])){
        DatasetSpecificScores[BroadClusterTypes[j]][[1]][k] <- Total_Mean
        if(is.na(Total_Mean)){DatasetSpecificScores[BroadClusterTypes[j]][[1]][k] <- 0}
      }
    }
  }
  ## The Final UCell score (TotalMinusEndo) for AD will be the upregulated scores minus the downregulated scores for all cell types except for endothelial cells.
  DatasetSpecificScores$TotalMinusEndo = DatasetSpecificScores$Astrocytes+DatasetSpecificScores$Oligodendrocytes+DatasetSpecificScores$Microglia+DatasetSpecificScores$Oligodendrocyte_Precursor_Cell+DatasetSpecificScores$Glut_ExcitatoryNeuron+DatasetSpecificScores$GABA_InhibitoryNeuron
  DatasetSpecificScores$Dataset = Datasets[i]
  write.table(DatasetSpecificScores, file=paste0("AD_UCellScores_UpRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testing_",Datasets[i],"_Dataset_Normalized2.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

###########################
### If there is no score present, then just set it to the mean of all other values (or to 0 if all of the data is not present)
for(i in 1:length(Datasets)){
  DatasetSpecificScores=read.table(paste0("AD_UCellScores_DownRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testing_",Datasets[i],"_Dataset_Normalized.txt"),header=T)
  for(j in 1:length(BroadClusterTypes)){
    Total_Mean = mean(DatasetSpecificScores[BroadClusterTypes[j]][[1]], na.rm=TRUE)
    for(k in 1:nrow(DatasetSpecificScores)){
      if(is.na(DatasetSpecificScores[BroadClusterTypes[j]][[1]][k])){
        DatasetSpecificScores[BroadClusterTypes[j]][[1]][k] <- Total_Mean
        if(is.na(Total_Mean)){DatasetSpecificScores[BroadClusterTypes[j]][[1]][k] <- 0}
      }
    }
  }
  ## The Final UCell score (TotalMinusEndo) for AD will be the upregulated scores minus the downregulated scores for all cell types except for endothelial cells.
  DatasetSpecificScores$TotalMinusEndo = DatasetSpecificScores$Astrocytes+DatasetSpecificScores$Oligodendrocytes+DatasetSpecificScores$Microglia+DatasetSpecificScores$Oligodendrocyte_Precursor_Cell+DatasetSpecificScores$Glut_ExcitatoryNeuron+DatasetSpecificScores$GABA_InhibitoryNeuron
  DatasetSpecificScores$Dataset = Datasets[i]
  write.table(DatasetSpecificScores, file=paste0("AD_UCellScores_DownRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testing_",Datasets[i],"_Dataset_Normalized2.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}










### This script gets the UCell scores for all datasets except a single dataset (the dataset held out is then tested based on the model created from all of the other datasets).
library(Seurat)
library(dplyr)
library(UCell)


#### Obtaining meta-data
Datasets = c("Mathys","Grubman","Lau","Morabito","Zhou","Leng_EC","Leng_SFG","OteroGarcia","Yang1","Yang2","Gerrits_OC","Gerrits_OTC","Smith_SSC","Smith_EC","Sadick","Gabitto_SeaAD","Barker","Sayed","Hoffman","Fujita")
BroadClusterTypes = c("Oligodendrocytes","Astrocytes","Oligodendrocyte_Precursor_Cell","Glut_ExcitatoryNeuron","Endothelial_Cell","GABA_InhibitoryNeuron","Microglia")
for(i in 1:length(Datasets)){
  PatientMetaData=read.table(paste(Datasets[i],"_Metadata_byIndividual.txt",sep=""),header=T)
  assign(paste("MetaData_",Datasets[i],sep=""),PatientMetaData)
}

#### Seurat pseudo-bulked files by average expression (not aggregate)
for(i in 1:length(Datasets)){
  currentTest = readRDS(paste0("avg_exp_",as.character(Datasets[i]),"_Mean_PredictedCellTypeL3.rds"))
  assign(paste("avg_exp_Mean_PredictedCellTypeL3_",Datasets[i],sep=""),currentTest)
}



for(l in 1:length(Datasets)){
  PredictedCellTypeL3_UpReg <- read.table(paste("SignificantGenes_Up.txt",sep=""),header=T)
  PredictedCellTypeL3_UpRegList = split(PredictedCellTypeL3_UpReg$Gene,PredictedCellTypeL3_UpReg$BroadClusterType)
  
  #### Goal: Use the genes to make UCell Scores for all other datasets.
  ### For each dataset have a table with individuals from all other datasets, 
  Datasets_Subset = Datasets[!(Datasets %in% Datasets[l])]
  FinalDatasetScoresTable = data.frame(1:1)
  FinalDatasetScoresTable[,2:(4+length(BroadClusterTypes))]=0
  colnames(FinalDatasetScoresTable)=c("Dataset","Individual","Disease_Severity","Diagnosis",BroadClusterTypes)
  for(z in 1:length(Datasets_Subset)){
    currentTest <- get(paste("avg_exp_Mean_PredictedCellTypeL3_",Datasets_Subset[z],sep=""))
    currentTest <- AddModuleScore_UCell(currentTest, features=PredictedCellTypeL3_UpRegList,maxRank=16000)
    UScoreMatrix = data.frame(1:length(unique(currentTest$Individual)))
    patients = unique(currentTest$Individual)
    UScoreMatrix[,1]= patients
    UScoreMatrix[,2:(length(BroadClusterTypes)+1)]=0
    colnames(UScoreMatrix)=c("Individual",BroadClusterTypes)
    for(j in 1:length(BroadClusterTypes)){
      if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets_Subset[z],][BroadClusterTypes[j]][[1]]==1){
        avg_temp <- subset(currentTest, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
        for(k in 1:length(patients)){
          avg_temp_Ind = avg_temp@meta.data[avg_temp$Individual==patients[k],]
          if(nrow(avg_temp_Ind)>0){
            CurrentName = paste0(BroadClusterTypes[j],"_UCell")
            UScoreMatrix[k,(j+1)]=avg_temp_Ind[CurrentName]
            ### If gene is not present, need to set UCell score to NA
            if(is.na(PredictedCellTypeL3_UpReg[PredictedCellTypeL3_UpReg$BroadClusterType==BroadClusterTypes[j],]$Gene[1])){UScoreMatrix[k,(j+1)]=NA}
          }
          if(nrow(avg_temp_Ind)==0){
            UScoreMatrix[k,(j+1)]=NA
          }
        }
      }
      if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets_Subset[z],][BroadClusterTypes[j]][[1]]==0){
        UScoreMatrix[,(j+1)]=NA
      }
    }
    PatientMetaData <- get(paste("MetaData_",Datasets_Subset[z],sep=""))
    UScoreMatrix_Expanded = cbind(UScoreMatrix,PatientMetaData)
    UScoreMatrix_Expanded$Dataset = Datasets_Subset[z]
    DatasetSpecificScores=UScoreMatrix_Expanded[c("Dataset","Individual","Disease_Severity","Diagnosis",BroadClusterTypes)]
    FinalDatasetScoresTable=rbind(FinalDatasetScoresTable,DatasetSpecificScores)
  }
  FinalDatasetScoresTable=FinalDatasetScoresTable[-1,]
  write.table(FinalDatasetScoresTable, file=paste0("AD_UCellScores_UpRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testingAllDatasets_LeavingOut_",Datasets[l],"_Dataset.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
  
  ####### Get DownRegulated Genes
  PredictedCellTypeL3_DownReg <- read.table(paste("SignificantGenes_Down.txt",sep=""),header=T)
  PredictedCellTypeL3_DownRegList = split(PredictedCellTypeL3_DownReg$Gene,PredictedCellTypeL3_DownReg$BroadClusterType)
  
  Datasets_Subset = Datasets[!(Datasets %in% Datasets[l])]
  FinalDatasetScoresTable = data.frame(1:1)
  FinalDatasetScoresTable[,2:(4+length(BroadClusterTypes))]=0
  colnames(FinalDatasetScoresTable)=c("Dataset","Individual","Disease_Severity","Diagnosis",BroadClusterTypes)
  for(z in 1:length(Datasets_Subset)){
    currentTest <- get(paste("avg_exp_Mean_PredictedCellTypeL3_",Datasets_Subset[z],sep=""))
    currentTest <- AddModuleScore_UCell(currentTest, features=PredictedCellTypeL3_DownRegList,maxRank=16000)
    UScoreMatrix = data.frame(1:length(unique(currentTest$Individual)))
    patients = unique(currentTest$Individual)
    UScoreMatrix[,1]= patients
    UScoreMatrix[,2:(length(BroadClusterTypes)+1)]=0
    colnames(UScoreMatrix)=c("Individual",BroadClusterTypes)
    for(j in 1:length(BroadClusterTypes)){
      if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets_Subset[z],][BroadClusterTypes[j]][[1]]==1){
        avg_temp <- subset(currentTest, subset=predicted.celltype.l3==as.character(BroadClusterTypes[j]))
        for(k in 1:length(patients)){
          avg_temp_Ind = avg_temp@meta.data[avg_temp$Individual==patients[k],]
          if(nrow(avg_temp_Ind)>0){
            CurrentName = paste0(BroadClusterTypes[j],"_UCell")
            UScoreMatrix[k,(j+1)]=avg_temp_Ind[CurrentName]
            ### If gene is not present, need to set UCell score to NA
            if(is.na(PredictedCellTypeL3_DownReg[PredictedCellTypeL3_DownReg$BroadClusterType==BroadClusterTypes[j],]$Gene[1])){UScoreMatrix[k,(j+1)]=NA}
          }
          if(nrow(avg_temp_Ind)==0){
            UScoreMatrix[k,(j+1)]=NA
          }
        }
      }
      if(PresenceofDataTable[PresenceofDataTable$Dataset==Datasets_Subset[z],][BroadClusterTypes[j]][[1]]==0){
        UScoreMatrix[,(j+1)]=NA
      }
    }
    PatientMetaData <- get(paste("MetaData_",Datasets_Subset[z],sep=""))
    UScoreMatrix_Expanded = cbind(UScoreMatrix,PatientMetaData)
    UScoreMatrix_Expanded$Dataset = Datasets_Subset[z]
    DatasetSpecificScores=UScoreMatrix_Expanded[c("Dataset","Individual","Disease_Severity","Diagnosis",BroadClusterTypes)]
    FinalDatasetScoresTable=rbind(FinalDatasetScoresTable,DatasetSpecificScores)
  }
  FinalDatasetScoresTable=FinalDatasetScoresTable[-1,]
  write.table(FinalDatasetScoresTable, file=paste0("AD_UCellScores_DownRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testingAllDatasets_LeavingOut_",Datasets[l],"_Dataset.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

### Normalize all scores.
##### UpReg
for(i in 1:length(Datasets)){
  DatasetSpecificScores=read.table(paste0("AD_UCellScores_UpRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testingAllDatasets_LeavingOut_",Datasets[i],"_Dataset.txt"),header=T)
  for(j in 1:length(BroadClusterTypes)){
    ClusterMean = mean(DatasetSpecificScores[BroadClusterTypes[j]][[1]], na.rm=TRUE)
    ClusterSum = sum(DatasetSpecificScores[BroadClusterTypes[j]][[1]], na.rm=TRUE)
    ## Only Normalize if there are some scores there for that cell type
    if(ClusterSum>0){
      ClusterMin = min(DatasetSpecificScores[BroadClusterTypes[j]][[1]],na.rm=TRUE)
      ClusterMax = max(DatasetSpecificScores[BroadClusterTypes[j]][[1]],na.rm=TRUE)
      ClusterRange = ClusterMax-ClusterMin
      DatasetSpecificScores[BroadClusterTypes[j]]=(DatasetSpecificScores[BroadClusterTypes[j]]-ClusterMin)/ClusterRange
    }
  }
  write.table(DatasetSpecificScores, file=paste0("AD_UCellScores_UpRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testingAllDatasets_LeavingOut_",Datasets[i],"_Dataset_Normalized.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}
##### DownReg
### Normalize all scores.
for(i in 1:length(Datasets)){
  DatasetSpecificScores=read.table(paste0("AD_UCellScores_DownRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testingAllDatasets_LeavingOut_",Datasets[i],"_Dataset.txt"),header=T)
  for(j in 1:length(BroadClusterTypes)){
    ClusterMean = mean(DatasetSpecificScores[BroadClusterTypes[j]][[1]], na.rm=TRUE)
    ClusterSum = sum(DatasetSpecificScores[BroadClusterTypes[j]][[1]], na.rm=TRUE)
    ## Only Normalize if there are some scores there for that cell type
    if(ClusterSum>0){
      ClusterMin = min(DatasetSpecificScores[BroadClusterTypes[j]][[1]],na.rm=TRUE)
      ClusterMax = max(DatasetSpecificScores[BroadClusterTypes[j]][[1]],na.rm=TRUE)
      ClusterRange = ClusterMax-ClusterMin
      DatasetSpecificScores[BroadClusterTypes[j]]=(DatasetSpecificScores[BroadClusterTypes[j]]-ClusterMin)/ClusterRange
    }
  }
  write.table(DatasetSpecificScores, file=paste0("AD_UCellScores_DownRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testingAllDatasets_LeavingOut_",Datasets[i],"_Dataset_Normalized.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}


###########################
### If there is no score present, then just set it to the mean of all other values (or to 0 if all of the data is not present)
###UpReg
for(i in 1:length(Datasets)){
  DatasetSpecificScores=read.table(paste0("AD_UCellScores_UpRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testingAllDatasets_LeavingOut_",Datasets[i],"_Dataset_Normalized.txt"),header=T)
  for(j in 1:length(BroadClusterTypes)){
    Total_Mean = mean(DatasetSpecificScores[BroadClusterTypes[j]][[1]], na.rm=TRUE)
    ## Loop through each individual: replace NA with the mean of that cell type's score.
    for(k in 1:nrow(DatasetSpecificScores)){
      if(is.na(DatasetSpecificScores[BroadClusterTypes[j]][[1]][k])){
        DatasetSpecificScores[BroadClusterTypes[j]][[1]][k] <- Total_Mean
        if(is.na(Total_Mean)){DatasetSpecificScores[BroadClusterTypes[j]][[1]][k] <- 0}
      }
    }
  }
  ## The Final UCell score (TotalMinusEndo) for AD will be the upregulated scores minus the downregulated scores for all cell types except for endothelial cells.
  DatasetSpecificScores$TotalMinusEndo = DatasetSpecificScores$Astrocytes+DatasetSpecificScores$Oligodendrocytes+DatasetSpecificScores$Microglia+DatasetSpecificScores$Oligodendrocyte_Precursor_Cell+DatasetSpecificScores$Glut_ExcitatoryNeuron+DatasetSpecificScores$GABA_InhibitoryNeuron
  write.table(DatasetSpecificScores, file=paste0("AD_UCellScores_UpRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testingAllDatasets_LeavingOut_",Datasets[i],"_Dataset_Normalized2.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

###########################
### If there is no score present, then just set it to the mean of all other values (or to 0 if all of the data is not present)
#### DownReg
for(i in 1:length(Datasets)){
  DatasetSpecificScores=read.table(paste0("AD_UCellScores_DownRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testingAllDatasets_LeavingOut_",Datasets[i],"_Dataset_Normalized.txt"),header=T)
  for(j in 1:length(BroadClusterTypes)){
    Total_Mean = mean(DatasetSpecificScores[BroadClusterTypes[j]][[1]], na.rm=TRUE)
    for(k in 1:nrow(DatasetSpecificScores)){
      if(is.na(DatasetSpecificScores[BroadClusterTypes[j]][[1]][k])){
        DatasetSpecificScores[BroadClusterTypes[j]][[1]][k] <- Total_Mean
        if(is.na(Total_Mean)){DatasetSpecificScores[BroadClusterTypes[j]][[1]][k] <- 0}
      }
    }
  }
  ## The Final UCell score (TotalMinusEndo) for AD will be the upregulated scores minus the downregulated scores for all cell types except for endothelial cells.
  DatasetSpecificScores$TotalMinusEndo = DatasetSpecificScores$Astrocytes+DatasetSpecificScores$Oligodendrocytes+DatasetSpecificScores$Microglia+DatasetSpecificScores$Oligodendrocyte_Precursor_Cell+DatasetSpecificScores$Glut_ExcitatoryNeuron+DatasetSpecificScores$GABA_InhibitoryNeuron
  write.table(DatasetSpecificScores, file=paste0("AD_UCellScores_DownRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testingAllDatasets_LeavingOut_",Datasets[i],"_Dataset_Normalized2.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}












library(pROC)

### Get Regression Scores
RegressionScores = data.frame(1:length(Datasets))
RegressionScores[,2:3]=0
colnames(RegressionScores)=c("Dataset","LogisticRegression_AUC","LinearRegression_CorrActvsPred")

for(i in 1:length(Datasets)){
  ### Get dataset to create the model
  DatasetSpecificScores_UpReg=read.table(paste0("AD_UCellScores_UpRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testingAllDatasets_LeavingOut_",Datasets[i],"_Dataset_Normalized2.txt"),header=T)
  DatasetSpecificScores_DownReg=read.table(paste0("AD_UCellScores_DownRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testingAllDatasets_LeavingOut_",Datasets[i],"_Dataset_Normalized2.txt"),header=T)
  RegressionScores[i,1]=Datasets[i]
  FinalDataset=DatasetSpecificScores_UpReg[c("Dataset","Individual","Disease_Severity","Diagnosis")]
  FinalDataset$Diagnosis[FinalDataset$Diagnosis=='AD'] <- '1'
  FinalDataset$Diagnosis[FinalDataset$Diagnosis=='Control'] <- '0'
  FinalDataset$Diagnosis[!(FinalDataset$Diagnosis==0 | FinalDataset$Diagnosis==1)] = NA
  FinalDataset$TotalMinusEndoTotal = DatasetSpecificScores_UpReg$TotalMinusEndo - DatasetSpecificScores_DownReg$TotalMinusEndo
  
  FinalDataset_forModel = FinalDataset
  FinalDataset_forModel2 = FinalDataset_forModel[!(is.na(FinalDataset_forModel$Diagnosis)),]
  FinalDataset_forModel2$Diagnosis = as.numeric(FinalDataset_forModel2$Diagnosis)
  
  FinalDataset_forModel=FinalDataset_forModel2[FinalDataset_forModel2$Diagnosis==1,]
  
  #### Get dataset to test:
  DatasetSpecificScores_UpReg=read.table(paste0("AD_UCellScores_UpRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testing_",Datasets[i],"_Dataset_Normalized2.txt"),header=T)
  DatasetSpecificScores_DownReg=read.table(paste0("AD_UCellScores_DownRegBySeverityScore_Top_",PercentageofTopDatasets,"_ofDatasets_testing_",Datasets[i],"_Dataset_Normalized2.txt"),header=T)
  RegressionScores[i,1]=Datasets[i]
  FinalDataset=DatasetSpecificScores_UpReg[c("Dataset","Individual","Disease_Severity","Diagnosis")]
  FinalDataset$Diagnosis[FinalDataset$Diagnosis=='AD'] <- '1'
  FinalDataset$Diagnosis[FinalDataset$Diagnosis=='Control'] <- '0'
  FinalDataset$Diagnosis[!(FinalDataset$Diagnosis==0 | FinalDataset$Diagnosis==1)] = NA
  FinalDataset$TotalMinusEndoTotal = DatasetSpecificScores_UpReg$TotalMinusEndo - DatasetSpecificScores_DownReg$TotalMinusEndo
  
  ### Test dataset
  FinalDataset2 = FinalDataset[!(is.na(FinalDataset$Diagnosis)),]
  FinalDataset2$Diagnosis = as.numeric(FinalDataset2$Diagnosis)
  
  FinalDataset=FinalDataset2[FinalDataset2$Diagnosis==1,]
  
  ### Make models
  ## Logistic Regression
  testLogReg <- glm(Diagnosis ~ TotalMinusEndoTotal,data=FinalDataset_forModel2, family="binomial")
  predicted <- testLogReg %>% predict(FinalDataset2, type = "response")
  testauc_log = auc(FinalDataset2$Diagnosis, predicted)
  
  RegressionScores[i,2]=testauc_log[1]
  
  ### Linear Regression (only testing cases)
  testLinReg <- lm(formula = Disease_Severity ~ TotalMinusEndoTotal, data= FinalDataset_forModel)
  predicted <- testLinReg %>% predict(FinalDataset, type = "response")
  testtable = data.frame(1:length(predicted))
  testtable[,1]=predicted
  testtable[,2]=FinalDataset$Disease_Severity
  RegressionScores[i,3]=cor(testtable[,1],testtable[,2])
}
write.table(RegressionScores,paste0("AD_RegressionScores_LeaveOneOutTestFullPredictionvsTestModel_Top_",PercentageofTopDatasets,"_ofDatasets.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE)

