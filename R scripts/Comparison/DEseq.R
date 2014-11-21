DESeq_pile <- function(count_table,normal,control) {
  library(DESeq)
  count_table<-count_table[,c(seq(1,normal),(dim(count_table)[2]/2)+seq(1,control))]
  condition<-c(substr(colnames(count_table),1,1))
  cds<-newCountDataSet(count_table,condition)
  cds<-estimateSizeFactors(cds)
  if(normal==1 &&control==1){
  cds<-estimateDispersions(cds,method='blind',sharingMode="fit-only")
  }else{
    cds<-estimateDispersions(cds)
  }
  output_DESeq<-nbinomTest(cds,levels(factor(condition))[1],levels(factor(condition))[2])
  return(output_DESeq)
}
