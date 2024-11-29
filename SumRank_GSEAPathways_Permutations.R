# DatasetNames is a vector of strings of dataset names. BroadClusterTypes is a vector of cell types. CommonGenes is a vector of genes held in common with all datasets.
# PresenceofDataTable is a table detailing whether a cell type has data for each cell type.
# This code requires the user to have files named DatasetName_CellType_DifferentialExpression_SuffixofDifferentialExpressionOutput.txt in the directory
# SuffixofDifferentialExpressionOutput is the end of the name you gave your differential expression file (e.g. "UpReg").
# This function will output multiple files including CellType_CombinedSignedNegLogPValranksNormalized_DirwinHallPValues_Top_ProportionofTopDatasetNames_ofDatasetNames_GSEAPathwayRanks.txt, which includes p-values for all pathways based on their reproducibility.
# m_df_gene2term is a data frame of 2 columns, the first of gene set names and the second of gene symbols for which pathway the genes are in.
SumRank_GSEAPathways_Permutations <- function(DatasetNames,BroadClusterTypes,SuffixofDifferentialExpressionOutput="",CommonGenes,ProportionofTopDatasetNames,PresenceofDataTable,directory,m_df_gene2term,FinalPathwayList_Table){
#### Step 1: Rank by Signed NegLog p-value
for(i in 1:length(DatasetNames)){
	for(j in 1:length(BroadClusterTypes)){
		if(PresenceofDataTable[PresenceofDataTable$Dataset==DatasetNames[i],][BroadClusterTypes[j]][[1]]==1){
			if(!(SuffixofDifferentialExpressionOutput=="")){currentTest <- read.table(paste(DatasetNames[i],"_",BroadClusterTypes[j],"_DifferentialExpression_",as.character(SuffixofDifferentialExpressionOutput),".txt",sep=""),header=T)}
			if((SuffixofDifferentialExpressionOutput=="")){currentTest <- read.table(paste(DatasetNames[i],"_",BroadClusterTypes[j],"_DifferentialExpression.txt",sep=""),header=T)}
			## CommonGenes is a vector of genes held in common with all datasets
			currentTest = currentTest[currentTest$Gene %in% CommonGenes,]
			currentTest$NegLogPVal = -log10(currentTest$p_val)
			currentTest$SignedNegLogPVal=currentTest$NegLogPVal
			## Set the SignedNegLogPVal to be negative if the log2FC is negative.
			currentTest$SignedNegLogPVal[currentTest$avg_log2FC<0] <- -currentTest$SignedNegLogPVal[currentTest$avg_log2FC<0]
			currentTest = currentTest[order(currentTest$SignedNegLogPVal,decreasing=TRUE),]
			currentTest$SignedNegLogPVal_rank = rank(-currentTest$SignedNegLogPVal)
      assign(paste0(DatasetNames[i],"_",BroadClusterTypes[j],"_SignedNegLogPVal_ranked"),currentTest)
      		}
    	}
  }
  
for(i in 1:length(DatasetNames)){
  for(j in 1:length(BroadClusterTypes)){
    if(PresenceofDataTable[PresenceofDataTable$Dataset==DatasetNames[i],][BroadClusterTypes[j]][[1]]==1){
      currentTest <- get(paste0(DatasetNames[i],"_",BroadClusterTypes[j],"_SignedNegLogPVal_ranked"))
      currentTest_geneList = setNames(currentTest$SignedNegLogPVal,currentTest$Gene)
      m_df_gene2term = m_df_gene2term[m_df_gene2term$gene_symbol %in% CommonGenes,]
      temporary = GSEA(geneList = currentTest_geneList, TERM2GENE=m_df_gene2term,pvalueCutoff=1)
      PathwaysDataframe = as.data.frame(temporary)
      PathwaysDataframe = data.frame(ID=PathwaysDataframe$ID,pvalue=PathwaysDataframe$pvalue,enrichmentScore=PathwaysDataframe$enrichmentScore)
      PathwaysDataframe$SignedNegLogPVal = -log10(PathwaysDataframe$pvalue)
      PathwaysDataframe$SignedNegLogPVal[PathwaysDataframe$enrichmentScore<0] <- -PathwaysDataframe$SignedNegLogPVal[PathwaysDataframe$enrichmentScore<0]
      PathwaysDataframe = PathwaysDataframe[order(PathwaysDataframe$SignedNegLogPVal,decreasing=TRUE),]
      PathwaysDataframe$SignedNegLogPVal_rank = 1:nrow(PathwaysDataframe)
      assign(paste0(DatasetNames[i],"_",BroadClusterTypes[j],"_SignedNegLogPVal_ranked_GSEAPathwayRanks"),PathwaysDataframe)
    }
  }
}

## Subset all datasets to FinalPathway list.
for(i in 1:length(DatasetNames)){
  for(j in 1:length(BroadClusterTypes)){
    if(PresenceofDataTable[PresenceofDataTable$Dataset==DatasetNames[i],][BroadClusterTypes[j]][[1]]==1){
      currentTest <- get(paste0(DatasetNames[i],"_",BroadClusterTypes[j],"_SignedNegLogPVal_ranked_GSEAPathwayRanks"))
      currentTest = currentTest[currentTest$ID %in% FinalPathwayList_Table[,1],]
      currentTest$SignedNegLogPVal_rank=1:nrow(currentTest)
      assign(paste0(DatasetNames[i],"_",BroadClusterTypes[j],"_SignedNegLogPVal_ranked_GSEAPathwayRanks_Subsetted"),currentTest)
}}}

#############
SumRankTable=data.frame(1:nrow(FinalPathwayList_Table))
SumRankTable[,1]=FinalPathwayList_Table[,1][order(as.character(FinalPathwayList_Table[,1]))]
SumRankTable[,2:(length(DatasetNames)+2)]=0
colnames(SumRankTable)=c("Gene",DatasetNames,"SumRankAcrossDatasetNames")

### Add all pvalranks together
## Loop through each cell type
for(j in 1:length(BroadClusterTypes)){
  ## Loop through each dataset
  for(i in 1:length(DatasetNames)){
    if(PresenceofDataTable[PresenceofDataTable$Dataset==DatasetNames[i],][BroadClusterTypes[j]][[1]]==1){
      currentTest <- get(paste0(DatasetNames[i],"_",BroadClusterTypes[j],"_SignedNegLogPVal_ranked_GSEAPathwayRanks_Subsetted"))
      currentTest = currentTest[order(currentTest$ID),]
      SumRankTable[,i+1]=(currentTest$SignedNegLogPVal_rank-1)/(nrow(FinalPathwayList_Table)-1)
    }
    if(PresenceofDataTable[PresenceofDataTable$Dataset==DatasetNames[i],][BroadClusterTypes[j]][[1]]==0){
      SumRankTable[,i+1]=NA
    }
  }
  SumRankTable[,(length(DatasetNames)+2)]=rowSums(SumRankTable[,2:(length(DatasetNames)+1)],na.rm=TRUE)
  assign(paste0(BroadClusterTypes[j],"_CombinedSignedNegLogPValranksNormalized_GSEAPathwayRanks"),SumRankTable)
}

###### Make a new table with the best (lowest) values for each pathway:
for(j in 1:length(BroadClusterTypes)){
  SumRankTable <- get(paste0(BroadClusterTypes[j],"_CombinedSignedNegLogPValranksNormalized_GSEAPathwayRanks"))
  ## Get a dataframe without the first and last columns to find the number of relevant DatasetNames
  test=as.numeric(SumRankTable[1,][,-1])
  test=head(test, -1)
  NumberofRelevantDatasetNames=length(na.exclude(test))
  NumberofRelevantDatasetNamesForRounding = NumberofRelevantDatasetNames
  ## Make Smaller tables that fit SixtyPercent of the number of relevant datasets.
  SumRankTable2=SumRankTable
  ## Find out how many NA columns there are:
  NumberNAColumns = sum(is.na(test))
  sorted_rankings = t(apply(SumRankTable[,2:(ncol(SumRankTable)-1)],MARGIN=1,sort))
  ## If there are NA columns, then add them to the end.
  if(NumberNAColumns>0){
    TabletoAdd = data.frame(1:nrow(sorted_rankings))
    TabletoAdd[,1:NumberNAColumns]=NA
    sorted_rankings=cbind(sorted_rankings,TabletoAdd)
  }
  ## Get the top half of datasets by sorting.
  SumRankTable2[,2:(ncol(SumRankTable2)-1)] = sorted_rankings
  Small_SumRankTable=SumRankTable2[,1:(NumberofRelevantDatasetNamesForRounding*as.numeric(ProportionofTopDatasets)+1)]
  Small_SumRankTable[,((NumberofRelevantDatasetNamesForRounding*as.numeric(ProportionofTopDatasets))+2)]=rowSums(Small_SumRankTable[,2:((NumberofRelevantDatasetNamesForRounding*as.numeric(ProportionofTopDatasets))+1)])
  assign(paste0(BroadClusterTypes[j],"_CombinedSignedNegLogPValranksNormalized_Top_",as.character(ProportionofTopDatasets),"_ofDatasets_GSEAPathwayRanks"),Small_SumRankTable)
}

###### Get p-values for each gene in each cell type by Irwin Hall distribution.
for(j in 1:length(BroadClusterTypes)){
	SumRankTable <- get(paste0(BroadClusterTypes[j],"_CombinedSignedNegLogPValranksNormalized_Top_",as.character(ProportionofTopDatasets),"_ofDatasets_GSEAPathwayRanks"))
  PVal_DirwinHallTable=data.frame(1:nrow(FinalPathwayList_Table))
  PVal_DirwinHallTable[,1]=FinalPathwayList_Table[,1][order(as.character(FinalPathwayList_Table[,1]))]
  #### If the sum rank is greater than half, you need to set it to half so it won't be considered significant.
  NumberofRelevantDatasetNames=ncol(SumRankTable)-2
  SumRankTable[,ncol(SumRankTable)][SumRankTable[,ncol(SumRankTable)] > (NumberofRelevantDatasetNames/2)] <- NumberofRelevantDatasetNames/2
  for(i in 1:nrow(FinalPathwayList_Table)){
    PVal_DirwinHallTable[i,2]=dirwin.hall(SumRankTable[i,ncol(SumRankTable)],(NumberofRelevantDatasetNames))
  }
  colnames(PVal_DirwinHallTable)=c("Gene","P_Val")
  PVal_DirwinHallTable$PlotPoint=1:nrow(PVal_DirwinHallTable)
  PVal_DirwinHallTable$NegLogPValue = -log10(PVal_DirwinHallTable$P_Val)
	if(!(SuffixofDifferentialExpressionOutput=="")){write.table(PVal_DirwinHallTable, file=paste(BroadClusterTypes[j],"_CombinedSignedNegLogPValranksNormalized_DirwinHallPValues_Top_",as.character(ProportionofTopDatasets),"_ofDatasets_",as.character(SuffixofDifferentialExpressionOutput),"_GSEAPathwayRanks.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)}
	if((SuffixofDifferentialExpressionOutput=="")){write.table(PVal_DirwinHallTable, file=paste(BroadClusterTypes[j],"_CombinedSignedNegLogPValranksNormalized_DirwinHallPValues_Top_",as.character(ProportionofTopDatasets),"_ofDatasets_GSEAPathwayRanks.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)}
}
}
