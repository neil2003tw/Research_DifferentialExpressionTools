DEcomparison<-function(count_table,AnswerDE,normalnum,controlnum){
  
  time_matrix<-matrix(nrow = 5,ncol = 1)
  rownames(time_matrix)<-c('DEGseq','DESeq','edgeR','limma','baySeq')
  colnames(time_matrix)<-'time'
  
  ##DEpipeline
  print(paste0('start rep count ',normalnum,'..'))
  
  print('DEGseq')
  time_matrix[1,1]<-system.time(output_DEGseq<-DEGseq_pipe(count_table,normalnum,controlnum))[3]
  sig_DEGseq<-output_DEGseq[ output_DEGseq$p.value<0.05 & !is.na(output_DEGseq$p.value) , ]
  
  print('DESeq')
  time_matrix[2,1]<-system.time(output_DESeq<-DESeq_pile(count_table,normalnum,controlnum))[3]
  output_DESeq<-output_DESeq[order(output_DESeq$pval,decreasing = F),]
  sig_DEseq<-output_DESeq[ output_DESeq$pval<0.05, ]
  
  print('edgeR')
  time_matrix[3,1]<-system.time(output_edgeR<-edgeR_pile(count_table,normalnum,controlnum))[3]
  output_edgeR<-output_edgeR$table[order(output_edgeR$table$PValue,decreasing = F),]
  sig_edgeR<-output_edgeR[ output_edgeR$PValue<0.05, ]
  
  print('limma')
  if (normalnum!=1 && controlnum!=1){
    time_matrix[4,1]<-system.time(output_limma<-limma_pile(count_table,normalnum,controlnum))[3]
    output_limma<-as.data.frame(output_limma$p.value[order(output_limma$p.value[,2],decreasing = F),])
    sig_limma<-output_limma[ output_limma$groupB<0.05, ]
  }else{
    print('condition unfit')
  }
  
  print('baySeq')
  if (normalnum!=1 && controlnum!=1){
    time_matrix[5,1]<-system.time(output_baySeq<-baySeq_pipe(count_table,normalnum,controlnum))[3]
    sig_baySeq<-output_baySeq[ output_baySeq$Likelihood>0.5, ]
  }
  
  
  ##Set summarise matrix
  
  generalmatrix<-matrix(nrow = 5,ncol = 4)
  colnames(generalmatrix)<-c('Sensitivity','Specificity','Precision','Accuracy')
  rownames(generalmatrix)<-c('DEGseq','DESeq','edgeR','limma','baySeq')
  
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
  
  
  if (normalnum!=1 && controlnum!=1){
    limmaTP<-sum(row.names(sig_limma) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=4 ])
    limmaFP<-sum(row.names(sig_limma) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==4 ])
    generalmatrix<-matrix_filledin(generalmatrix,4,limmaTP,limmaFP,CP,CN)
  }else{
    print('limma no s/s values')
  }
  
  bayTP<-sum(row.names(sig_baySeq) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)!=4 ])
  bayFP<-sum(row.names(sig_baySeq) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==4 ])
  generalmatrix<-matrix_filledin(generalmatrix,5,bayTP,bayFP,CP,CN)
  
  ##check find out ratio of different type of DE
  
  DE_based_matrix<-matrix(nrow = 5,ncol = 3)
  colnames(DE_based_matrix)<-c('highDE','midDE','lowDE')
  rownames(DE_based_matrix)<-c('DEGseq','DESeq','edgeR','limma','baySeq')
  DEcase<-c('highDE','midDE','lowDE')
  
  for(i in seq(1,3)){
    DE_based_matrix[1,i]<-sum(sig_DEGseq$GeneNames %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
    DE_based_matrix[2,i]<-sum(sig_DEseq$id %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
    DE_based_matrix[3,i]<-sum(row.names(sig_edgeR) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
    if (normalnum!=1 && controlnum!=1){
      DE_based_matrix[4,i]<-sum(row.names(sig_limma) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
    }else{
      print('limma gone again')
    }
    DE_based_matrix[5,i]<-sum(row.names(sig_baySeq) %in% AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])/length(AnswerDE$gene_name[ as.numeric(AnswerDE$DE)==i ])
  }
  
  
  Alltestingdata<-cbind(generalmatrix,DE_based_matrix,time_matrix)
  
  
  draw_all_ven<-function(){  
    library(grid)
    DE_DEGseq<-merge(data.frame(gene_name=sig_DEGseq$GeneNames),AnswerDE,by='gene_name')
    DE_DESeq<-merge(data.frame(gene_name=sig_DEseq$id),AnswerDE,by='gene_name')
    DE_edgeR<-merge(data.frame(gene_name=row.names(sig_edgeR)),AnswerDE,by='gene_name')
    DE_limma<-merge(data.frame(gene_name=row.names(sig_limma)),AnswerDE,by='gene_name')
    if (!file.exists('img')){
      dir.create('img')
    } 
    wid=350
    hei=300
    grid.newpage()
    cQ<-quad_Ven(DE_DEGseq[DE_DEGseq$DE!='unDE',],DE_DESeq[DE_DESeq$DE!='unDE',],DE_edgeR[DE_edgeR$DE!='unDE',],DE_limma[DE_limma$DE!='unDE',])
    png(filename = paste0('img/',normalnum,controlnum,'CQ.png'),width = wid,height = hei)
    grid.draw(cQ)
    dev.off()
    grid.newpage()
    hQ<-quad_Ven(DE_DEGseq[DE_DEGseq$DE=='highDE',],DE_DESeq[DE_DESeq$DE=='highDE',],DE_edgeR[DE_edgeR$DE=='highDE',],DE_limma[DE_limma$DE=='highDE',])
    png(filename = paste0('img/',normalnum,controlnum,'HQ.png'),width = wid,height = hei)
    grid.draw(hQ)
    dev.off()
    grid.newpage()
    mQ<-quad_Ven(DE_DEGseq[DE_DEGseq$DE=='midDE',],DE_DESeq[DE_DESeq$DE=='midDE',],DE_edgeR[DE_edgeR$DE=='midDE',],DE_limma[DE_limma$DE=='midDE',])
    png(filename = paste0('img/',normalnum,controlnum,'MQ.png'),width = wid,height = hei)
    grid.draw(mQ)
    dev.off()
    grid.newpage()
    lQ<-quad_Ven(DE_DEGseq[DE_DEGseq$DE=='lowDE',],DE_DESeq[DE_DESeq$DE=='lowDE',],DE_edgeR[DE_edgeR$DE=='lowDE',],DE_limma[DE_limma$DE=='lowDE',])
    png(filename = paste0('img/',normalnum,controlnum,'LQ.png'),width = wid,height = hei)
    grid.draw(lQ)
    dev.off()
  }
  if (normalnum!=1 && controlnum!=1){
    draw_all_ven()
  }
  return(Alltestingdata)
}

