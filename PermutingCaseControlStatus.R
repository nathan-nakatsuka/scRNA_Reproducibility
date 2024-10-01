# This code permutes cases and controls. It requires the following: DatasetName (a string with the name of the Dataset in PresenceofDataTable)
# PresenceofDataTable (a table indicating whether the dataset has cells for each cell type), CellTypeName (a string indicating what is the cell type label, such as "predicted.celltype.l1")
# BroadClusterTypes (a vector of strings naming all cell types), CaseName (a string naming what cases are called in the dataset, such as "COVID-19"), and ControlName (a string naming what controls are called in the dataset, such as "healthy").
# Note: for this code you must call the individual ID as "patient" and the disease name as "disease_status_standard".
PermuteCaseControl <- function(DatasetName,avg_exp_Dataset,PresenceofDataTable,CellTypeName,BroadClusterTypes,CaseName,ControlName){
	individuals = unique(avg_exp_Dataset$patient)
        individuals = na.omit(individuals)
        Individualtodisease_status_standardTable = data.frame(1:length(individuals))
        Individualtodisease_status_standardTable[,1]=individuals
        Individualtodisease_status_standardTable[,2:3]=0
        ### Find their normal disease_status_standard and Sex
        for(r in 1:nrow(Individualtodisease_status_standardTable)){
          Individualtodisease_status_standardTable[r,2]=avg_exp_Dataset@meta.data[avg_exp_Dataset@meta.data$patient==Individualtodisease_status_standardTable[r,1],]$disease_status_standard[1]
        }
        ###Permuting the Individuals' case/control status
        d=0
        while(d==0){
          rm(.Random.seed)
          Individualtodisease_status_standardTable[,3]=sample(Individualtodisease_status_standardTable[,2],length(Individualtodisease_status_standardTable[,2]),replace=FALSE)
          rm(.Random.seed)
          avg_exp_Dataset$disease_status_standard_Permuted = avg_exp_Dataset$disease_status_standard
          ### Put new permuted disease_status_standard labels to object
          for(s in 1:nrow(avg_exp_Dataset@meta.data)){
            avg_exp_Dataset$disease_status_standard_Permuted[s]=Individualtodisease_status_standardTable[match(avg_exp_Dataset@meta.data$patient[s],Individualtodisease_status_standardTable[,1]),3]
          }
          ### Test if all cell types that are present have at least 2 of each 
          Case_metadata = avg_exp_Dataset@meta.data[avg_exp_Dataset@meta.data$disease_status_standard_Permuted==CaseName,]
          Control_metadata = avg_exp_Dataset@meta.data[avg_exp_Dataset@meta.data$disease_status_standard_Permuted==ControlName,]
          TestClusters={}
          	for(j in 1:length(BroadClusterTypes)){
          		if(PresenceofDataTable[PresenceofDataTable$Dataset==DatasetName,(j+1)]==1){TestClusters=append(TestClusters,BroadClusterTypes[j])}
          	}
          CaseTest = data.frame(1:length(TestClusters))
          ControlTest = data.frame(1:length(TestClusters))
          for(w in 1:length(TestClusters)){
            CaseTest[w,1]=nrow(Case_metadata[Case_metadata[[CellTypeName]]==TestClusters[w],])
            ControlTest[w,1]=nrow(Control_metadata[Control_metadata[[CellTypeName]]==TestClusters[w],])
          }
          if(all(CaseTest[,1]>1) & all(ControlTest[,1]>1)){d=1}
        }
	return(avg_exp_Dataset)
	}
