source('DEGseq.R')
source('DEseq.R')
source('limma.R')
source('baySeq.R')
source('edgeR.R')
source('counttable_merge.R')

count_table<-counttable_merged('../../Normal distributed case/')
AnswerDE<-read.csv('../../Normal distributed case/allDElist.csv')
savetohere<-function(data,filename){
  path<-paste0('../../Normal distributed case/',filename)
  write.table(data,path)
}

##DEpipeline
output_DESeq<-DESeq_pile(count_table,normal = 5,control = 5)
output_limma<-limma_pile(count_table,normal = 5,control = 5)
output_edgeR<-edgeR_pile(count_table,normal = 5,control = 5)
output_baySeq<-baySeq_pipe(count_table,normal = 5,control = 5)
output_DEGseq<-DEGseq_pipe(count_table,normal = 5,control = 5)


##data cleaning
output_DESeq<-output_DESeq[order(output_DESeq$pval,decreasing = F),]
output_edgeR<-output_edgeR$table[order(output_edgeR$table$PValue,decreasing = F),]
output_limma<-as.data.frame(output_limma$p.value[order(output_limma$p.value[,2],decreasing = F),])
  #output_baySeq

##filter out p>0.05
sig_DEGseq<-output_DEGseq[ output_DEGseq$p.value<0.05 & !is.na(output_DEGseq$p.value) , ]
sig_DEseq<-output_DESeq[ output_DESeq$pval<0.05, ]
sig_edgeR<-output_edgeR[ output_edgeR$PValue<0.05, ]
sig_limma<-output_limma[ output_limma$groupB<0.05, ]


##add group

##adjust answer matrix
AnswerDE<-AnswerDE[!is.na(AnswerDE$DE),]
levels(AnswerDE$DE)<-c('highDE','highDE','lowDE','lowDE','midDE','midDE','unDE')



##Set summarise matrix

generalmatrix<-matrix(nrow = 4,ncol = 4)
rownames(generalmatrix)<-c('Sensitivity','Specificity','Precision','Accuracy')
colnames(generalmatrix)<-c('DEGseq','DESeq','edgeR','limma')

matrix_filledin<-function(matrix,row,TP,FP,CP=CP,CN=CN){
  matrix[1,row]<-TP/CP
  matrix[2,row]<-(CN-FP)/CN
  matrix[3,row]<-TP/(TP+FP)
  matrix[4,row]<-(TP+(CN-FP))/(CP+CN)
  return(matrix)
}

CP<-length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=4 ] )
CN<-length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==4 ] )

DEGTP<-sum(sig_DEGseq$GeneNames %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=4 ])
DEGFP<-sum(sig_DEGseq$GeneNames %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==4 ])
generalmatrix<-matrix_filledin(generalmatrix,1,DEGTP,DEGFP,CP,CN)


DETP<-sum(sig_DEseq$id %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=4 ])
DEFP<-sum(sig_DEseq$id %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==4 ])
generalmatrix<-matrix_filledin(generalmatrix,2,DETP,DEFP,CP,CN)

edgeTP<-sum(row.names(sig_edgeR) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=4 ])
edgeFP<-sum(row.names(sig_edgeR) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==4 ])
generalmatrix<-matrix_filledin(generalmatrix,3,edgeTP,edgeFP,CP,CN)

limmaTP<-sum(row.names(sig_limma) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=4 ])
limmaFP<-sum(row.names(sig_limma) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==4 ])
generalmatrix<-matrix_filledin(generalmatrix,4,limmaTP,limmaFP,CP,CN)

savetohere(generalmatrix,'summary_all.txt')
rm(list=ls(pattern = '[PN]$'))



##check find out ratio of different type of DE

DE_based_matrix<-matrix(nrow = 3,ncol = 4)
rownames(DE_based_matrix)<-c('highDE','midDE','lowDE')
colnames(DE_based_matrix)<-c('DEGseq','DESeq','edgeR','limma')
DEcase<-c('highDE','midDE','lowDE')

for(i in seq(1,3)){
  DE_based_matrix[i,1]<-sum(sig_DEGseq$GeneNames %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
  DE_based_matrix[i,2]<-sum(sig_DEseq$id %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
  DE_based_matrix[i,3]<-sum(row.names(sig_edgeR) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
  DE_based_matrix[i,4]<-sum(row.names(sig_limma) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
}

savetohere(DE_based_matrix,'summary_DEbase.txt')


##HighDE summarise
highDEmatrix<-matrix(nrow = 4,ncol = 4)
rownames(highDEmatrix)<-c('Sensitivity','Specificity','Precision','Accuracy')
colnames(highDEmatrix)<-c('DEGseq','DESeq','edgeR','limma')


