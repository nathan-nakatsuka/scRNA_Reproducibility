# Produces a Manhattan plot from calibrated p-values ordered by genes alphabetically.
# CalibratedPValuesTable is a table of genes and their calibrated p-values comparing the real data p-values with p-values obtained from Permutations that can be obtained from the CalibratePValueswithPermutations function.
# OtherNegLogPValueCutoff is a cutoff that can be used that is usually lower than the significant value threshold (this can be obtained by cross-validation to see what maximized reproducibility).
# TopValueCutoff is the top negative log10 p-value allowed (values of infinite, which is obtained if the value is above that of every permutation value, will be set to this for plotting purposes).
# Desiredggtitle is the title of the ggplot. jitter_amount is the amount of jitter to add to points that are at the top due to having infinite negative log10 p-values.
MakeManhattanPlot <- function(CalibratedPValuesTable, OtherNegLogPValueCutoff,TopValueCutoff, Desiredggtitle,jitter_amount=0.01){
  CalibratedPValuesTable$PlotPoint = 1:nrow(CalibratedPValuesTable)
  CalibratedPValuesTable$New_NegLogPValue = -log10(CalibratedPValuesTable$PVal)
  # Generate jittered values for each infinite entry
  jittered_values <- as.numeric(TopValueCutoff) + jitter(rep(0, sum(is.infinite(CalibratedPValuesTable$New_NegLogPValue))), amount = jitter_amount)
  # Replace infinite values with the jittered values
  CalibratedPValuesTable[is.infinite(CalibratedPValuesTable$New_NegLogPValue), 5] <- jittered_values
  SignificantGenes = CalibratedPValuesTable[CalibratedPValuesTable$PVal_BH<0.05,]
  SignificantGenes2 = CalibratedPValuesTable[CalibratedPValuesTable$New_NegLogPValue>as.numeric(OtherNegLogPValueCutoff),]
  SignificantGenes2=SignificantGenes2[!(SignificantGenes2$Gene %in% SignificantGenes$Gene),]
  bad = CalibratedPValuesTable[CalibratedPValuesTable$PVal_BH>0.05,]
  Badcutoff = max(bad$New_NegLogPValue)+0.00001
  plot(CalibratedPValuesTable$PlotPoint,CalibratedPValuesTable$New_NegLogPValue,ylim=c(0,(as.numeric(TopValueCutoff)+1)),pch=20,xaxt='n',xlab="Genes",ylab=expression("-log"[10]*"(p-value)"),main=Desiredggtitle,cex.main=0.6)
  abline(h=Badcutoff, lty=3,col="red",lwd=1.5)
  if(nrow(SignificantGenes)>0){text(x=SignificantGenes$PlotPoint,y=SignificantGenes$New_NegLogPValue,labels=SignificantGenes$Gene,col="red",cex=0.5)}
  if(Badcutoff>as.numeric(OtherNegLogPValueCutoff)){abline(h=as.numeric(OtherNegLogPValueCutoff), lty=3,col="orange",lwd=1.5)}
  if(nrow(SignificantGenes2)>0){text(x=SignificantGenes2$PlotPoint,y=SignificantGenes2$New_NegLogPValue,labels=SignificantGenes2$Gene,col="orange",cex=0.5)}
}
