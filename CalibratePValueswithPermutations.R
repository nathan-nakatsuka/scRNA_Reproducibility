# Code like this can be used to combine the permutation results.
# Combine the permutation results
TotalNumberofPermutations = 100
ComparisonTable = data.frame(1:TotalNumberofPermutations)
for(z in 1:TotalNumberofPermutations){
    setwd(paste0("/home/mydirectory/Permutation",as.character(PermutationNumber)))
    PVal_DirwinHallTable <- read.table(paste(ClusterofInterest,"_CombinedSignedNegLogPValranksNormalized_DirwinHallPVals_Top_",as.character(ProportionofDatasetstoUse),"_ofDatasets_Permutation",as.character(z),".txt",sep=""),header=T)
    ComparisonTable[((z-1)*length(CommonGenes)+1):(z*length(CommonGenes)),1] = PVal_DirwinHallTable$NegLogPValue
}
write.table(ComparisonTable,"/home/mydirectory/PermutationComparisonTable.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)



CalibratePValueswithPermutations <- function(NumberofPermutations,ComparisonTable){
	PVal_DirwinHallTable <- read.table(paste(BroadClusterTypes[j],"_COVID_DESeq2_Pseudobulk_All_ControlvsCOVID_DownReg_CombinedSignedNegLogPValranksNormalized_DirwinHallPVals_Top",args[1],"Percentof16Datasets_AllData.txt",sep=""),header=T)
	FinalTable = data.frame(1:length(CommonGenes))
	FinalTable[,1]=CommonGenes
	FinalTable[,2:3]=0
	ecdf_fun = ecdf(ComparisonTable[,1])
	FinalTable[,2]=1-sapply(X = PVal_DirwinHallTable[,4], FUN = ecdf_fun)
	FinalTable[,3]=p.adjust(FinalTable[,2],method="BH")
	colnames(FinalTable)=c("Gene","PVal","PVal_BH")
	write.table(FinalTable, file=paste("/gpfs/commons/home/nnakatsuka/COVID/PermutationTests/Mean_SDTables/",BroadClusterTypes[j],"_COVID_DESeq2_Pseudobulk_All_ControlvsCOVID_DownReg_CombinedSignedNegLogPValranksNormalized_DirwinHallPVals_Top",args[1],"Percentof16Datasets_AllData_PermutationThresholdPValues.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

