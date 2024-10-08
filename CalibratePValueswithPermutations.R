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


# This function returns a table of genes and their calibrated p-values comparing the real data p-values with p-values obtained from Permutations.
# PVal_DirwinHallTable is the output from the SumRank function.
# ComparisonTable is a table listing all negativelogpvalues from SumRank done on permutations
# CommonGenes is a vector of genes held in common with all datasets that can be obtained with the GetCommonGenes function.
CalibratePValueswithPermutations <- function(CommonGenes,ComparisonTable,PVal_DirwinHallTable){
	FinalTable = data.frame(1:length(CommonGenes))
	FinalTable[,1]=CommonGenes
	FinalTable[,2:3]=0
	ecdf_fun = ecdf(ComparisonTable[,1])
	FinalTable[,2]=1-sapply(X = PVal_DirwinHallTable$NegLogPValue, FUN = ecdf_fun)
	FinalTable[,3]=p.adjust(FinalTable[,2],method="BH")
	colnames(FinalTable)=c("Gene","PVal","PVal_BH")
	return(FinalTable)
}

