
PermuteSex <- function(DatasetNames, PresenceofDataTable,CellTypeName,BroadClusterTypes){
  for(j in 1:length(BroadClusterTypes)){
    for(i in 1:length(DatasetNames)){
      currentTest <- get(paste("avg_exp_DESeq2_PredictedCellTypeL1_2_",DatasetNames[i],sep=""))
      currentTest <- subset(currentTest, disease_status_standard %in% c("COVID-19","healthy"))
		  currentTest[["RNA"]]@counts <- as.matrix(currentTest[["RNA"]]@counts)+1
      currentTest[["RNA"]]@counts <- round(currentTest[["RNA"]]@counts)
        if(PresenceofDataTable[PresenceofDataTable$Dataset==DatasetNames[i],][BroadClusterTypes[j]][[1]]==1){
		      currentTest <- subset(currentTest, subset=CellTypeName==as.character(BroadClusterTypes[j]))
      individuals = unique(currentTest$patient)
      individuals = na.omit(individuals)
      IndividualtoSexTable = data.frame(1:length(individuals))
      IndividualtoSexTable[,1]=individuals
      IndividualtoSexTable[,2:3]=0
      ### Find their normal disease_status_standard and Sex
      for(r in 1:nrow(IndividualtoSexTable)){
        IndividualtoSexTable[r,2]=currentTest@meta.data[currentTest@meta.data$patient==IndividualtoSexTable[r,1],]$disease_status_standard[1]
        IndividualtoSexTable[r,3]=currentTest@meta.data[currentTest@meta.data$patient==IndividualtoSexTable[r,1],]$sex_standard2[1]
      }
      ### Randomly permute the sexes within controls using the same set of sexes (so you end up with the same number of each sex within the controls)
        rm(.Random.seed)
		IndividualtoSexTable[IndividualtoSexTable[,2]=="healthy",3]=sample(IndividualtoSexTable[IndividualtoSexTable[,2]=="healthy",3],length(IndividualtoSexTable[IndividualtoSexTable[,2]=="healthy",3]),replace=FALSE)
        rm(.Random.seed)
        IndividualtoSexTable[IndividualtoSexTable[,2]=="COVID-19",3]=sample(IndividualtoSexTable[IndividualtoSexTable[,2]=="COVID-19",3],length(IndividualtoSexTable[IndividualtoSexTable[,2]=="COVID-19",3]),replace=FALSE)
        currentTest$sex_standard2_Permuted = currentTest$sex_standard2
        ### Put new permuted sex labels to object
        for(s in 1:nrow(currentTest@meta.data)){
          currentTest$sex_standard2_Permuted[s]=IndividualtoSexTable[match(currentTest@meta.data$patient[s],IndividualtoSexTable[,1]),3]
        }
}
          
