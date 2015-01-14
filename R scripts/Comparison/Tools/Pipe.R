DEGseq_pipe<-function(count_table,normal,control){
  library(DEGseq)
  count_table<-count_table[,c(seq(1,normal),(dim(count_table)[2]/2)+seq(1,control))]
  write.table(data.frame(row.names(count_table),count_table),'count_table.txt',quote = F,sep='\t',row.names=F,col.names = F)
  geneMatrix1<-readGeneExp(file='count_table.txt', geneCol=1, valCol=c(1+seq(1,normal)),header = F,sep = '\t')
  geneMatrix2<-readGeneExp(file='count_table.txt', geneCol=1, valCol=c(1+normal+seq(1,control)),header = F,sep = '\t')
  DEGexp(geneExpMatrix1=geneMatrix1, geneCol1=1, expCol1=c(1+seq(1,normal)), groupLabel1="A",geneExpMatrix2=geneMatrix2, geneCol2=1, expCol2=c(1+seq(1,control)), groupLabel2="B",method="MARS",outputDir = '.',pValue = 0.05)
  system('rm -r output/ output.html count_table.txt')
  DEGseq_out<-read.table('output_score.txt',header = T)
  system('rm -r output_score.txt')
  return(DEGseq_out)
}

####
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

###
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

###
limma_pile<-function(count_table,normal,control){
  library(edgeR)
  library(limma)
  count_table<-count_table[,c(seq(1,normal),(dim(count_table)[2]/2)+seq(1,control))]
  condition<-c(substr(colnames(count_table),1,1))
  dge<-DGEList(count_table,group = factor(condition)) 
  dge<-calcNormFactors(dge)
  v<-voom(dge,plot=F)
  fit<-lmFit(v)
  fit<-eBayes(fit)
  return(fit)
}

###
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