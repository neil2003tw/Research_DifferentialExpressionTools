baySeq_pipe<-function(count_table,normal,control){
library(baySeq)
library(snow)
cl<-makeCluster(4,'SOCK')
count_table<-count_table[,c(seq(1,normal),(dim(count_table)[2]/2)+seq(1,control))]
CD <- new("countData", data = as.matrix(count_table), replicates = colnames(count_table), groups = list(NDE=c(rep(1,normal+control)),DE=c(rep(1,normal),rep(2,control))))
libsizes(CD) <- getLibsizes(CD)
CD <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", cl = cl)
CD <- getLikelihoods.NB(CD, pET = 'BIC', cl = cl)
return(topCounts(CD,group=2,number = 1000000))
}