limma_pile<-function(count_table,normal,control){
library(edgeR)
library(limma)
count_table<-count_table[,c(seq(1,normal),(dim(count_table)[2]/2)+seq(1,control))]
condition<-c(substr(colnames(count_table),1,1))
dge<-DGEList(count_table,group = factor(condition)) 
dge<-calcNormFactors(dge)
v<-voom(dge,plot=T)
fit<-lmFit(v)
fit<-eBayes(fit)
return(fit)
}