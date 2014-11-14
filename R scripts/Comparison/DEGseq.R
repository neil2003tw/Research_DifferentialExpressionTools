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