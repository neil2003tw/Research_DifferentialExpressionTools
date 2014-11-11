limma_pile<-function(count_table){
library(edgeR)
library(limma)
dge<-DGEList(count_table,group = factor(condition)) 
dge<-calcNormFactors(dge)
v<-voom(dge,plot=TRUE)
fit<-lmFit(v)
fit<-eBayes(fit)
return(fit)
}