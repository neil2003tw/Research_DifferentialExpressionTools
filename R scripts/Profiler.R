library(caret)

flux_profile_maker<-function(name,profiledata,replicate=1){
  
  makeprofile<-function(x,group='A'){
    if(group=='A'){
      gene_name[x]<-0
      gene_name[[x]][highADE]<-250*10
      gene_name[[x]][highBDE]<-250*1
      gene_name[[x]][midADE]<-250*5
      gene_name[[x]][midBDE]<-250*1
      gene_name[[x]][lowADE]<-250*2
      gene_name[[x]][lowBDE]<-250*1
      gene_name[[x]][unDE]<-250*1
      
    }else if(group=='B'){
      gene_name[x]<-0
      gene_name[[x]][highBDE]<-250*10
      gene_name[[x]][highADE]<-250*1
      gene_name[[x]][midBDE]<-250*5
      gene_name[[x]][midADE]<-250*1
      gene_name[[x]][lowBDE]<-250*2
      gene_name[[x]][lowADE]<-250*1
      gene_name[[x]][unDE]<-250*1
    }
    merged_profile<-merge(flux_profile,gene_name,by.x = 'chr_position',by.y = 'gene_name')
    merged_profile[[x]]<-sapply(merged_profile[[x]],function(x) as.integer(x+rnorm(1,100,50)))
    merged_profile$ratio<-0
    merged_profile$ratio<-sapply(merged_profile[[x]], function(y) y/sum(merged_profile[[x]],na.rm = TRUE))
    file<-paste0(name,'/profiled_',x,'.pro')
    merged_profile<-data.frame(merged_profile[,1],merged_profile[,2],merged_profile[,3],merged_profile[,4],merged_profile[,6],merged_profile[,5])
    write.table(merged_profile,file=file,sep = '\t',row.names=FALSE,quote=FALSE)
    print(x)
  }
  
  if (!file.exists(name)){
    dir.create(name)
  } 
  
  t<-read.table(profiledata,sep='\t')
  flux_profile<-data.frame(t$V1,t$V2,t$V3,t$V4)
  colnames(flux_profile)<-c('chr_position','gene_id','region','gene_length')
  gene_name<-levels(flux_profile$chr_position)

  set.seed(11111)
  index_gene<-seq(1:length(gene_name))
  highDE<-sample(index_gene,1000)
  middleDE<-sample(index_gene[-highDE],2000)
  DEed<-c(highDE,middleDE)
  lowDE<-sample(index_gene[-DEed],3000)
  DEed<-c(DEed,lowDE)
  unDE<-sample(index_gene[-DEed],5000)


  HA<-createDataPartition(highDE,p=0.5)
  highADE<-highDE[HA$Resample1]
  highBDE<-highDE[-HA$Resample1]
  MA<-createDataPartition(middleDE,p=0.5)
  midADE<-middleDE[MA$Resample1]
  midBDE<-middleDE[-MA$Resample1]
  LA<-createDataPartition(lowDE,p=0.5)
  lowADE<-lowDE[LA$Resample1]
  lowBDE<-lowDE[-LA$Resample1]

  gene_name<-data.frame(gene_name)

  Agroup<-paste0(rep('A',replicate),seq(1:replicate))
  Bgroup<-paste0(rep('A',replicate),seq(1:replicate))
  print('READIED')
  for(i in 1:replicate){
   makeprofile(Agroup[i],'A')
  }
  for(i in 1:replicate){
   makeprofile(Bgroup[i],'B')
  }

  gene_name$DE<-NA
  gene_name$DE[highADE]<-'highADE'
  gene_name$DE[highBDE]<-'highBDE'
  gene_name$DE[midADE]<-'midADE'
  gene_name$DE[midBDE]<-'midBDE'
  gene_name$DE[lowADE]<-'lowADE'
  gene_name$DE[lowBDE]<-'lowBDE'
  gene_name$DE[unDE]<-'unDE'
  gene_name$DE<-factor(gene_name$DE)
  write.csv(gene_name,file=paste0(name,'/allDElist.csv'),row.names=FALSE)
}