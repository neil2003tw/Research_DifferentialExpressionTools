DEGseq_pipe<-function(count_table){
  library(DEGseq)
  write.table(data.frame(row.names(count_table),count_table),'count_table.txt',quote = F,sep='\t',row.names=F,col.names = F)
  geneMatrix1<-readGeneExp(file='count_table.txt', geneCol=1, valCol=c(2,3,4,5),header = F,sep = '\t')
  geneMatrix2<-readGeneExp(file='count_table.txt', geneCol=1, valCol=c(6,7,8,9),header = F,sep = '\t')
  DEGexp(geneExpMatrix1=geneMatrix1, geneCol1=1, expCol1=c(2,3,4,5), groupLabel1="A",geneExpMatrix2=geneMatrix2, geneCol2=1, expCol2=c(2,3,4,5), groupLabel2="B",method="MARS",outputDir = '.',pValue = 0.05)
  system('rm -r output/ output.html count_table.txt')
  DEGseq_out<-read.table('output_score.txt',header = T)
  return(DEGseq_out)
  system('rm -r output_score.txt')
}