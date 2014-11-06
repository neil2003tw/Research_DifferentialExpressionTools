DESeq_pile <- function(file, ) {
count_table<-read.csv(file ,row.names=1)
condition<-c(colnames(count_table))
cds<-newCountDataSet(count_table,condition)
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds,method='blind',sharingMode='fit-only')
}
