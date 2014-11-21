edgeR_pile<-function(count_table,normal,control){
library(edgeR)
count_table<-count_table[,c(seq(1,normal),(dim(count_table)[2]/2)+seq(1,control))]
condition<-c(substr(colnames(count_table),1,1))
y<-DGEList(count_table,group = factor(condition)) 
if(!(normal==1 &&control==1)){
  y<-calcNormFactors(y)
  y<-estimateCommonDisp(y)
  y<-estimateTagwiseDisp(y)
  output_edgeR<-exactTest(y)
}else{
  output_edgeR<-exactTest(y,dispersion = 0.4^2)
}
return(output_edgeR)
}
