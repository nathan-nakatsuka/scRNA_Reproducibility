
PermuteCaseControl <- function(DatasetNames,PresenceofDataTable,CellTypeName){
	for(i in 1:length(DatasetNames)){
		currentTest <- get(paste("avg_exp_DESeq2_PredictedCellTypeL1_2_",DatasetNames[i],sep=""))
        individuals = unique(currentTest$patient)
        individuals = na.omit(individuals)
        Individualtodisease_status_standardTable = data.frame(1:length(individuals))
        Individualtodisease_status_standardTable[,1]=individuals
        Individualtodisease_status_standardTable[,2:3]=0
        ### Find their normal disease_status_standard and Sex
        for(r in 1:nrow(Individualtodisease_status_standardTable)){
          Individualtodisease_status_standardTable[r,2]=currentTest@meta.data[currentTest@meta.data$patient==Individualtodisease_status_standardTable[r,1],]$disease_status_standard[1]
        }
        ###Permuting the Individuals' case/control status
        d=0
        while(d==0){
          rm(.Random.seed)
          Individualtodisease_status_standardTable[,3]=sample(Individualtodisease_status_standardTable[,2],length(Individualtodisease_status_standardTable[,2]),replace=FALSE)
          rm(.Random.seed)
          currentTest$disease_status_standard_Permuted = currentTest$disease_status_standard
          ### Put new permuted disease_status_standard labels to object
          for(s in 1:nrow(currentTest@meta.data)){
            currentTest$disease_status_standard_Permuted[s]=Individualtodisease_status_standardTable[match(currentTest@meta.data$patient[s],Individualtodisease_status_standardTable[,1]),3]
          }
          ### Test if all cell types that are present have at least 2 of each 
          AD_metadata = currentTest@meta.data[currentTest@meta.data$disease_status_standard_Permuted=="COVID-19",]
          Control_metadata = currentTest@meta.data[currentTest@meta.data$disease_status_standard_Permuted=="healthy",]
          TestClusters={}
          	for(j in 1:length(BroadClusterTypes)){
          		if(PresenceofDataTable[PresenceofDataTable$Dataset==DatasetNames[i],(j+1)]==1){TestClusters=append(TestClusters,BroadClusterTypes[j])}
          	}
          ADTest = data.frame(1:length(TestClusters))
          ControlTest = data.frame(1:length(TestClusters))
          for(w in 1:length(TestClusters)){
            ADTest[w,1]=nrow(AD_metadata[AD_metadata$predicted.celltype.l1_2==TestClusters[w],])
            ControlTest[w,1]=nrow(Control_metadata[Control_metadata$predicted.celltype.l1_2==TestClusters[w],])
          }
          if(all(ADTest[,1]>1) & all(ControlTest[,1]>1)){d=1}
        }
}
}
