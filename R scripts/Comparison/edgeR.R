edgeR_pile<-function(count_table){
library(edgeR)
condition<-c(substr(colnames(count_table),1,1))
y<-DGEList(count_table,group = factor(condition)) 
y<-calcNormFactors(y)
y<-estimateCommonDisp(y)
y<-estimateTagwiseDisp(y)
output_edgeR<-exactTest(y)
return(output_edgeR)
}
