limma_pile<-function(count_table){
library(edgeR)
library(limma)
condition<-c(substr(colnames(count_table),1,1))
dge<-DGEList(count_table,group = factor(condition)) 
dge<-calcNormFactors(dge)
v<-voom(dge,plot=T)
fit<-lmFit(v)
fit<-eBayes(fit)
return(fit)
}