# This function produces a table with the number of cases and controls and total number of individuals in each dataset.
# The code requires Seurat objects named DatasetName_Seurat, which have IDs labeled "patient" and disease status labeled "disease_status_standard".
# Datasets is a vector of dataset names. CaseName is a string naming what the cases are called (e.g. "COVID-19"). ControlName is a string naming what the controls are called (e.g. "healthy")
MakeDatasetNumberTable <- function(Datasets,CaseName,ControlName){
FinalTable = data.frame(1:length(Datasets))
FinalTable[,2:4]=0
colnames(FinalTable)=c("Dataset","NumCases","NumControls","TotalNum")
for(i in 1:length(Datasets)){
  temp = get(paste0(Datasets[i],"_Seurat"))
  temp_metadata = temp@meta.data[!duplicated(temp@meta.data$patient), ]
  FinalTable[i,1]=Datasets[i]
  FinalTable[i,2]=nrow(temp_metadata[temp_metadata$disease_status_standard==CaseName,])
  FinalTable[i,3]=nrow(temp_metadata[temp_metadata$disease_status_standard==ControlName,])
}
FinalTable[,4]=FinalTable[,2]+FinalTable[,3]
return(FinalTable)
}
