


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
output_DESeq<-DESeq_pile(count_table,normal = 1,control = 1)
output_limma<-limma_pile(count_table,normal = 1,control = 1)
output_edgeR<-edgeR_pile(count_table,normal = 1,control = 1)
output_baySeq<-baySeq_pipe(count_table,normal = 1,control = 1)
output_DEGseq<-DEGseq_pipe(count_table,normal = 1,control = 1)


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
colnames(generalmatrix)<-c('Sensitivity','Specificity','Precision','Accuracy')
rownames(generalmatrix)<-c('DEGseq','DESeq','edgeR','limma')

matrix_filledin<-function(matrix,row,TP,FP,CP=CP,CN=CN){
  matrix[row,1]<-TP/CP
  matrix[row,2]<-(CN-FP)/CN
  matrix[row,3]<-TP/(TP+FP)
  matrix[row,4]<-(TP+(CN-FP))/(CP+CN)
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

DE_based_matrix<-matrix(nrow = 4,ncol = 3)
rownames(DE_based_matrix)<-c('highDE','midDE','lowDE')
colnames(DE_based_matrix)<-c('DEGseq','DESeq','edgeR','limma')
DEcase<-c('highDE','midDE','lowDE')

for(i in seq(1,3)){
  DE_based_matrix[1,i]<-sum(sig_DEGseq$GeneNames %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
  DE_based_matrix[2,i]<-sum(sig_DEseq$id %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
  DE_based_matrix[3,i]<-sum(row.names(sig_edgeR) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
  DE_based_matrix[4,i]<-sum(row.names(sig_limma) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
}

savetohere(DE_based_matrix,'summary_DEbase.txt')


##Ven diagram
quad_Ven<-function(d1,d2,d3,d4){

library(VennDiagram)

cn1<-length(d1[ ,1 ])
cn2<-length(d2[ ,1 ])
cn3<-length(d3[ ,1 ])
cn4<-length(d4[ ,1 ])

cn12<-sum(d1[,1] %in% d2[,1])
cn13<-sum(d1[,1] %in% d3[,1])
cn14<-sum(d1[,1] %in% d4[,1])
cn23<-sum(d2[,1] %in% d3[,1])
cn24<-sum(d2[,1] %in% d4[,1])
cn34<-sum(d3[,1] %in% d4[,1])

cn123<-sum( d1[d1[,1] %in% d2[,1],1] %in% d3[,1])
cn124<-sum( d1[d1[,1] %in% d2[,1],1] %in% d4[,1])
cn134<-sum( d1[d1[,1] %in% d3[,1],1] %in% d4[,1])
cn234<-sum( d2[d2[,1] %in% d3[,1],1] %in% d4[,1])

cn1234<-sum( d1[d1[,1] %in% d2[,1],1] %in% d3[d3[,1] %in% d4[,1],1])

venn.plot<-draw.quad.venn(cn1,cn2,cn3,cn4,cn12,cn13,cn14,cn23,cn24,cn34,cn123,cn124,cn134,cn234,cn1234,
                          category=c('DEGseq','DESeq','edgeR','limma'),fill = c("orange", "red", "green", "blue"),
                          lty = "dashed",cex = 2, cat.cex = 1.5,cat.col = c("orange", "red", "green", "blue"))

}


DE_DEGseq<-merge(data.frame(gene_name=sig_DEGseq$GeneNames),AnswerDE,by='gene_name')
DE_DESeq<-merge(data.frame(gene_name=sig_DEseq$id),AnswerDE,by='gene_name')
DE_edgeR<-merge(data.frame(gene_name=row.names(sig_edgeR)),AnswerDE,by='gene_name')
DE_limma<-merge(data.frame(gene_name=row.names(sig_limma)),AnswerDE,by='gene_name')


wid=350
hei=300
cQ<-quad_Ven(DE_DEGseq[DE_DEGseq$DE!='unDE',],DE_DESeq[DE_DESeq$DE!='unDE',],DE_edgeR[DE_edgeR$DE!='unDE',],DE_limma[DE_limma$DE!='unDE',])
png(filename = 'CQ.png',width = wid,height = hei)
grid.draw(cQ)
dev.off()
grid.newpage()
hQ<-quad_Ven(DE_DEGseq[DE_DEGseq$DE=='highDE',],DE_DESeq[DE_DESeq$DE=='highDE',],DE_edgeR[DE_edgeR$DE=='highDE',],DE_limma[DE_limma$DE=='highDE',])
png(filename = 'HQ.png',width = wid,height = hei)
grid.draw(hQ)
dev.off()
grid.newpage()
mQ<-quad_Ven(DE_DEGseq[DE_DEGseq$DE=='midDE',],DE_DESeq[DE_DESeq$DE=='midDE',],DE_edgeR[DE_edgeR$DE=='midDE',],DE_limma[DE_limma$DE=='midDE',])
png(filename = 'MQ.png',width = wid,height = hei)
grid.draw(mQ)
dev.off()
grid.newpage()
lQ<-quad_Ven(DE_DEGseq[DE_DEGseq$DE=='lowDE',],DE_DESeq[DE_DESeq$DE=='lowDE',],DE_edgeR[DE_edgeR$DE=='lowDE',],DE_limma[DE_limma$DE=='lowDE',])
png(filename = 'LQ.png',width = wid,height = hei)
grid.draw(lQ)
dev.off()

