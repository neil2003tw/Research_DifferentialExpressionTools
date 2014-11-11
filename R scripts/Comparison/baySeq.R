baySeq_pipe<-function(count_table){
library(baySeq)
library(snow)
cl<-makeCluster(4,'SOCK')
replicates=length(colnames(count_table))/2
CD <- new("countData", data = as.matrix(count_table), replicates = colnames(count_table), groups = list(NDE=c(rep(1,2*replicates)),DE=c(rep(1,replicates),rep(2,replicates))))
libsizes(CD) <- getLibsizes(CD)
CD <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", cl = cl)
CD <- getLikelihoods.NB(CD, pET = 'BIC', cl = cl)
return(topCounts(CD,group=2,number = 1000000))
}