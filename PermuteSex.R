# This code permutes sex. It requires the following: DatasetName (a string with the name of the Dataset in PresenceofDataTable)
# PresenceofDataTable (a table indicating whether the dataset has cells for each cell type), CellTypeName (a string indicating what is the cell type label, such as "predicted.celltype.l1"),
# avg_exp_Dataset is a Seurat object that the user will use for differential expression. This is usually a pseudobulked object (or not pseudobulked if using a mixed model for differential expression).
# BroadClusterTypes (a vector of strings naming all cell types), CaseName (a string naming what cases are called in the dataset, such as "COVID-19"), and ControlName (a string naming what controls are called in the dataset, such as "healthy").
# Note: for this code you must call the sex as "sex_standard", the individual ID as "patient", and the disease name as "disease_status_standard" in the avg_exp_Dataset Seurat Object.
# This will output a new Seurat object with the column name sex_standard_Permuted

PermuteSex <- function(DatasetName,avg_exp_Dataset,PresenceofDataTable,CellTypeName,BroadClusterTypes,CaseName,ControlName){
	currentMetaData = avg_exp_Dataset@meta.data
	currentMetaData$Index = paste0(currentMetaData$patient,"_",currentMetaData[[CellTypeName]])
	currentMetaData$sex_standard_Permuted = NA
	for(j in 1:length(BroadClusterTypes)){
	        if(PresenceofDataTable[PresenceofDataTable$Dataset==DatasetName[i],][BroadClusterTypes[j]][[1]]==1){
			currentTest = currentMetaData[currentMetaData[[CellTypeName]]==as.character(BroadClusterTypes[j]),]
      			individuals = unique(currentTest$patient)
      			individuals = na.omit(individuals)
      			IndividualtoSexTable = data.frame(1:length(individuals))
      			IndividualtoSexTable[,1]=individuals
      			IndividualtoSexTable[,2:3]=0
      			### Find their normal cell type, disease_status_standard and Sex
      			for(r in 1:nrow(IndividualtoSexTable)){
        			IndividualtoSexTable[r,2]=currentTest[currentTest$patient==IndividualtoSexTable[r,1],]$disease_status_standard[1]
        			IndividualtoSexTable[r,3]=currentTest[currentTest$patient==IndividualtoSexTable[r,1],]$sex_standard[1]
      			}
      			### Randomly permute the sexes within controls using the same set of sexes (so you end up with the same number of each sex within the controls)
        		rm(.Random.seed)
			IndividualtoSexTable[IndividualtoSexTable[,2]==ControlName,3]=sample(IndividualtoSexTable[IndividualtoSexTable[,2]==ControlName,3],length(IndividualtoSexTable[IndividualtoSexTable[,2]==ControlName,3]),replace=FALSE)
        		rm(.Random.seed)
        		IndividualtoSexTable[IndividualtoSexTable[,2]==CaseName,3]=sample(IndividualtoSexTable[IndividualtoSexTable[,2]==CaseName,3],length(IndividualtoSexTable[IndividualtoSexTable[,2]==CaseName,3]),replace=FALSE)
        		currentTest$sex_standard_Permuted = currentTest$sex_standard
        		### Put new permuted sex labels to object
        		for(s in 1:nrow(currentTest)){
          			currentTest$sex_standard_Permuted[s]=IndividualtoSexTable[match(currentTest$patient[s],IndividualtoSexTable[,1]),3]
        			}
		}
		# Put the new permuted sex labels to the metadata
		currentMetaData[match(currentTest$Index,currentMetaData$Index)]$sex_standard_Permuted = currentTest[match(currentMetaData$Index,currentTest$Index)]$sex_standard_Permuted
	}
	# Put the metadata back in the original object
	avg_exp_Dataset@meta.data=currentMetaData
	return(avg_exp_Dataset)
}
