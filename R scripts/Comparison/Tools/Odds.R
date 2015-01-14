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


###
savetohere<-function(data,filename){
  path<-paste0('../../Normal distributed case/',filename)
  write.table(data,path)
}


###

counttable_merged<-function(directory){
  
  origindir<-getwd()  
  setwd(directory)

  data_name<-dir()
  data_name<-data_name[grepl('profiled',data_name)]
  condition<-substr(data_name,10,11)
  datas<-paste0('profile',condition)
  for(x in seq(1,length(condition))){
    assign(datas[x],read.csv(data_name[x],sep = ' ',col.names=c('Gene',condition[x])))
  }
  list_profile<-lapply(ls(pattern = 'profile[AB][0-9]'), function(x) get(x))
  count_table<-Reduce(function(x, y) merge(x, y, all=TRUE), list_profile)
  count_table[is.na(count_table)]<-0
  count_table<-data.frame(count_table[,-1],row.names=count_table[,1])
  rm(list=datas)
  
  setwd(origindir)
  
  return(count_table)
}

