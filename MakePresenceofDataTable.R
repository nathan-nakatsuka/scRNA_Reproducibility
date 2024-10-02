# This is code to make the PresenceofDataTable file, which indicates what cell types are present
# This code requires the presence of Seurat objects (pseudobulked or not) that are named avg_exp_DatasetName, where DatasetName is the name of the dataset.
# BroadClusterTypes is a vector of strings naming all cell types
# CellTypeLevel indicates the cell type level resolution (e.g. "predicted.celltype.l1")
MakePresenceofDataTable <- function(DatasetNames,BroadClusterTypes,CellTypeLevel){
  PresenceofDataTable = data.frame(1:length(DatasetNames))
  PresenceofDataTable[,1]=DatasetNames
  PresenceofDataTable[,2:(length(BroadClusterTypes)+1)]=0
  colnames(PresenceofDataTable)=c("Dataset",BroadClusterTypes)
  for(i in 1:length(DatasetNames)){
    currentTest <- get(paste("avg_exp_",DatasetNames[i],sep=""))
    for(j in 1:length(BroadClusterTypes)){
      odd=currentTest@meta.data[currentTest@meta.data[[CellTypeLevel]]==BroadClusterTypes[j],]
      if(nrow(odd)>0){PresenceofDataTable[i,(j+1)]=1}
      if(nrow(odd)==0){PresenceofDataTable[i,(j+1)]=0}
    }
  }
  return(PresenceofDataTable)
}
