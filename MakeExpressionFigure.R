# This function makes an expression figure showing the relative expression and standard deviation of cases and controls.
# Datasets is a vector of dataset names. CellTypeLevel indicates the cell type level resolution (e.g. "predicted.celltype.l1"). 
# CaseName is a string naming what cases are called in the dataset, such as "COVID-19". ControlName is a string naming what controls are called in the dataset, such as "healthy".
# GeneofInterest is a string of what gene you want to test. CellType is a string of the cell type (e.g. "CD4_T"). 
# DatasetNumbersTable is a table with the numbers of cases and controls that can be made with MakeDatasetNumberTable.
MakeExpressionFigure <- function(Datasets,CellTypeLevel,CaseName,ControlName,GeneofInterest,CellType,DatasetNumbersTable){
FinalTable = data.frame(1:(length(Datasets)))
FinalTable[,1]=Datasets
FinalTable[,2:5]=0
colnames(FinalTable)=c("Dataset","Case_MeanExpression","Case_SDExpression","Control_MeanExpression","Control_SDExpression")
for(i in 1:length(Datasets)){
  currentTest <- get(paste0("avg_exp_Mean_",Datasets[i]))
  gene1<- FetchData(currentTest, vars = GeneofInterest)
  gene1$disease_status_standard = currentTest$disease_status_standard
  gene1$CellType = currentTest@meta.data[[CellTypeLevel]]
  celltype_expression = gene1[gene1$CellType==ClusterofInterest,]
  celltype_expression_cases = celltype_expression[celltype_expression$disease_status_standard==CaseName,]
  celltype_expression_controls = celltype_expression[celltype_expression$disease_status_standard==ControlName,]
  FinalTable[i,2]=mean(celltype_expression_cases[,1],na.rm=T)
  FinalTable[i,3]=sd(celltype_expression_cases[,1],na.rm=T)
  FinalTable[i,4]=mean(celltype_expression_controls[,1],na.rm=T)
  FinalTable[i,5]=sd(celltype_expression_controls[,1],na.rm=T)
}

DatasetNumbersTable=DatasetNumbersTable[order(DatasetNumbersTable$Dataset),]
FinalTable=FinalTable[order(FinalTable$Dataset),]
FinalTable$TotalNum = DatasetNumbersTable$TotalNum
FinalTable2=FinalTable
temp1=min(FinalTable$Case_MeanExpression-FinalTable$Case_SDExpression)
temp2=max(FinalTable$Case_MeanExpression+FinalTable$Case_SDExpression)
temp3=min(FinalTable$Control_MeanExpression-FinalTable$Control_SDExpression)
temp4=max(FinalTable$Control_MeanExpression+FinalTable$Control_SDExpression)
temp5=min(temp1,temp3)
temp6=max(temp2,temp4)
p <- ggplot(data=FinalTable, aes(y=Case_MeanExpression,x=Control_MeanExpression))+geom_point(aes(size=TotalNum))+theme_classic()+theme(plot.title=element_text(size=20),axis.title=element_text(size=16),legend.title=element_text(size=16))+geom_abline(intercept=0,slope=1)+labs(x="Normalized Mean Expression in Controls",y="Normalized Mean Expresion in Cases")+ggtitle(paste0("Expression of ",GeneofInterest," in ",ClusterofInterest," \n in Cases and Controls"))+theme(plot.title = element_text(hjust=0.5))+xlim((temp5-0.0001),(temp6+0.0001))+ylim((temp5-0.0001),(temp6+0.00001))   
p+guides(size=guide_legend("Dataset Sample \n Size"))+geom_pointrange(aes(ymin=Case_MeanExpression-Case_SDExpression,ymax=Case_MeanExpression+Case_SDExpression), linewidth=.01)+geom_errorbarh(aes(xmin=Control_MeanExpression-Control_SDExpression, xmax=Control_MeanExpression+Control_SDExpression),height=0,linewidth=0.01)
}
