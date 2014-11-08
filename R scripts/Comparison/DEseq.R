DESeq_pile <- function(file, ) {
condition<-c(substr(colnames(count_table),1,1))
cds<-newCountDataSet(count_table,condition)
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)
output_DESeq<-nbinomTest(cds,levels(factor(condition))[1],levels(factor(condition))[2])
}
