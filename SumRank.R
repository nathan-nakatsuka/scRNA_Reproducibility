library(unifed)

GetCommonGenes <- function(DatasetNames,BroadClusterTypes,PresenceofDataTable){
  #Find a cell type with no missing data.
  CellTypeIndexwithNoMissing = which(colSums(PresenceofDataTable[2:ncol(PresenceofDataTable)])==nrow(PresenceofDataTable))[[1]]
  ## Get list of common genes 
  objs <- list()
  for(i in 1:length(DatasetNames)){
    currentTest <- read.table(paste(DatasetNames[i],"_",BroadClusterTypes[CellTypeIndexwithNoMissing],"_DifferentialExpression.txt",sep=""),header=T)
    objs[[i]] <- currentTest$Gene
  }
  CommonGenes = Reduce(intersect, objs)
  CommonGenes=CommonGenes[order(as.character(CommonGenes))]
  return(CommonGenes)
}

SumRank <- function(DatasetNames,BroadClusterTypes,CommonGenes,ProportionofTopDatasets,PresenceofDataTable,directory){
#### Step 1: Rank by Signed NegLog p-value
for(i in 1:length(DatasetNames)){
	for(j in 1:length(BroadClusterTypes)){
		if(PresenceofDataTable[PresenceofDataTable$Dataset==DatasetNames[i],][BroadClusterTypes[j]][[1]]==1){
			currentTest <- read.table(paste(DatasetNames[i],"_",BroadClusterTypes[j],"_DifferentialExpression.txt",sep=""),header=T)
			## CommonGenes is a vector of genes held in common with all datasets
			currentTest = currentTest[currentTest$Gene %in% CommonGenes,]
			currentTest$NegLogPVal = -log10(currentTest$p_val)
			currentTest$SignedNegLogPVal=currentTest$NegLogPVal
			## Set the SignedNegLogPVal to be negative if the log2FC is negative.
			currentTest$SignedNegLogPVal[currentTest$avg_log2FC<0] <- -currentTest$SignedNegLogPVal[currentTest$avg_log2FC<0]
			currentTest = currentTest[order(currentTest$SignedNegLogPVal,decreasing=TRUE),]
			currentTest$SignedNegLogPVal_rank = rank(-currentTest$SignedNegLogPVal)
			write.table(currentTest, file=paste(directory,DatasetNames[i],"_",BroadClusterTypes[j],"_SignedNegLogPVal_ranked.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
      }
    }
  }
	
#### Step 2: Add the signed neglog10p-value ranks of all datasets together
SumRankTable=data.frame(1:length(CommonGenes))
SumRankTable[,1]=CommonGenes[order(as.character(CommonGenes))]
## Loop through each cell type
for(j in 1:length(BroadClusterTypes)){
	for(i in 1:length(DatasetNames)){
		if(PresenceofDataTable[PresenceofDataTable$Dataset==DatasetNames[i],][BroadClusterTypes[j]][[1]]==1){
			currentTest <- read.table(paste(DatasetNames[i],"_",BroadClusterTypes[j],"_SignedNegLogPVal_ranked.txt",sep=""),header=T)
			## Order all datasets by Gene (so they can be matched properly)
			currentTest = currentTest[order(currentTest$Gene),]
			### Normalize the ranks
			SumRankTable[,i+1]=(currentTest$SignedNegLogPVal_rank-1)/(length(CommonGenes)-1)
		}
		if(PresenceofDataTable[PresenceofDataTable$Dataset==DatasetNames[i],][BroadClusterTypes[j]][[1]]==0){
			SumRankTable[,i+1]=NA
		}
	}
	SumRankTable[,(length(DatasetNames)+2)]=rowSums(SumRankTable[,2:(length(DatasetNames)+1)])
	colnames(SumRankTable)=c("Gene",DatasetNames[1:length(DatasetNames)],"SumRankAcrossDatasets")
	write.table(SumRankTable, file=paste(BroadClusterTypes[j],"_CombinedSignedNegLogPValranksNormalized.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

#### Step 3:  Take the top datasets (choose a percentage) for each gene.
for(j in 1:length(BroadClusterTypes)){
	SumRankTable <- read.table(paste(BroadClusterTypes[j],"_CombinedSignedNegLogPValranksNormalized.txt",sep=""),header=T)
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
	Small_SumRankTable=SumRankTable2[,1:(NumberofRelevantDatasetsForRounding*as.numeric(ProportionofTopDatasets)+1)]
	Small_SumRankTable[,((NumberofRelevantDatasetsForRounding*as.numeric(ProportionofTopDatasets))+2)]=rowSums(Small_SumRankTable[,2:((NumberofRelevantDatasetsForRounding*as.numeric(ProportionofTopDatasets))+1)])
	write.table(Small_SumRankTable, file=paste(BroadClusterTypes[j],"_CombinedSignedNegLogPValranksNormalized_Top_",as.character(ProportionofTopDatasets),"_ofDatasets.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}

#### Step 4: Get p-values for each gene by Irwin Hall distribution
for(j in 1:length(BroadClusterTypes)){
	SumRankTable <- read.table(paste(BroadClusterTypes[j],"_CombinedSignedNegLogPValranksNormalized_Top_",as.character(ProportionofTopDatasets),"_ofDatasets.txt",sep=""),header=F)
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
	write.table(PVal_DirwinHallTable, file=paste(BroadClusterTypes[j],"_CombinedSignedNegLogPValranksNormalized_DirwinHallPValues_Top_",as.character(ProportionofTopDatasets),"_ofDatasets.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}
}
