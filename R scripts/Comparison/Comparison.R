source('DEGseq.R')
source('DEseq.R')
source('limma.R')
source('baySeq.R')
source('edgeR.R')

##DEpipeline
output_DESeq<-DESeq_pile(count_table)
output_limma<-limma_pile(count_table)
output_edgeR<-edgeR_pile(count_table)
output_baySeq<-baySeq_pipe(count_table)
output_DEGseq<-DEGseq_pipe(count_table)


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

##read Golden standard
AnswerDE<-read.csv('../../Random/allDElist.csv')
AnswerDE<-AnswerDE[!is.na(AnswerDE$DE),]

CP<-length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=7 ])
CN<-length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==7] )


##Set final matrix
Finalmatrix<-matrix(nrow = 4,ncol = 4)
rownames(Finalmatrix)<-c('Sensitivity','Specificity','Precision','Accuracy')
colnames(Finalmatrix)<-c('DEGseq','DESeq','edgeR','limma')

DEGTP<-sum(sig_DEGseq$GeneNames %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=7 ])
DEGTN<-sum(sig_DEGseq$GeneNames %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==7 ])

DETP<-sum(sig_DEseq$id %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=7 ])
DETN<-sum(sig_DEseq$id %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==7 ])

edgeTP<-sum(row.names(sig_edgeR) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=7 ])
edgeTN<-sum(row.names(sig_edgeR) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==7 ])

limmaTP<-sum(row.names(sig_limma) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=7 ])
limmaTN<-sum(row.names(sig_limma) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=7 ])


