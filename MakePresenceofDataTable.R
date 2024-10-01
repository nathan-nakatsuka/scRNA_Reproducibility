# This is code to make the PresenceofDataTable file

MakePresenceofDataTable <- function(DatasetNames,BroadClusterTypes){
PresenceofDataTable = data.frame(1:length(DatasetNames))
PresenceofDataTable[,1]=DatasetNames
PresenceofDataTable[,2:(length(BroadClusterTypes)+1)]=0
colnames(PresenceofDataTable)=c("Dataset",BroadClusterTypes)

setwd("/brahms/nakatsukan/COVID/DiffExpGenes/DESeq2_Full_Test/PredictedCellTypesL1")
for(i in 1:length(Datasets)){
  for(j in 1:length(BroadClusterTypes)){
    currentTest <- get(paste("avg_exp_",Datasets[i],sep=""))
    SubsettedtoCellType = currentTest@meta.data[currentTest@meta.data$predicted.celltype.l1_2==BroadClusterTypes[j],]
    CasesSubset = SubsettedtoCellType[SubsettedtoCellType$disease_status_standard=="COVID-19",]
    ControlsSubset = SubsettedtoCellType[SubsettedtoCellType$disease_status_standard=="healthy",]
    if(nrow(CasesSubset)>1 & nrow(ControlsSubset)>1){
      avg_temp <- subset(currentTest, subset=predicted.celltype.l1_2==as.character(BroadClusterTypes[j]))
      avg_temp[["RNA"]]@counts = round(avg_temp[["RNA"]]@counts)
      Idents(avg_temp) <- "disease_status_standard"
      temp10 <- FindMarkers(avg_temp, ident.1="COVID-19",ident.2="healthy",test.use="DESeq2",group.by = 'disease_status_standard', logfc.threshold = 0, min.pct=0,min.cells.group=1)
      temp10$Gene = rownames(temp10)
      temp10 = temp10[!grepl('MT-',temp10$Gene),]
      write.table(temp10, file=paste(Datasets[i],"_",BroadClusterTypes[j],"_COVID_DESeq2_Pseudobulk_All_UpReg.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
      PresenceofDataTable[i,(j+1)]=1
    }
  }}
