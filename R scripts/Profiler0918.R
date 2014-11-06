library(caret)

flux_profile_maker<-function(name,profiledata,replicate=1){
  
  makeprofile<-function(){
      
      group_A<-paste0(rep('A',replicate),seq(1,replicate))
      group_B<-paste0(rep('B',replicate),seq(1,replicate))
      
 
      for(i in seq(1,replicate)){
        gene_name[group_A[i]]<-0
      }
      for(i in seq(1,replicate)){
        gene_name[group_B[i]]<-0
      }
      
      for(i in highADE){
        maxt<-rnorm(1,mean=2.331907,sd = 1.728949)
        while(maxt>7.251596){
          maxt<-rnorm(1,mean=2.331907,sd = 1.728949)
        }
        gene_name[i,seq(2,replicate+1)]<-maxt
      }
      for(i in highBDE){
        maxt<-rnorm(1,mean=2.331907,sd = 1.728949)
        while(maxt>7.251596){
          maxt<-rnorm(1,mean=2.331907,sd = 1.728949)
        }
        gene_name[i,seq(replicate+2,(2*replicate)+1)]<-maxt
      }
      
      
      for(i in midADE){
        maxt<-rnorm(1,mean=2.285371,sd = 1.697909)
        while(maxt>4.800113){
          maxt<-rnorm(1,mean=2.285371,sd = 1.697909)
        }
        gene_name[i,seq(2,replicate+1)]<-maxt
      }
      for(i in midBDE){
        maxt<-rnorm(1,mean=2.285371,sd = 1.697909)
        while(maxt>4.800113){
          maxt<-rnorm(1,mean=2.285371,sd = 1.697909)
        }
        gene_name[i,seq(replicate+2,(2*replicate)+1)]<-maxt
      }
      
      
      for(i in lowADE){
        maxt<-rnorm(1,mean=2.050947,sd = 1.498367)
        while(maxt>4.536703){
          maxt<-rnorm(1,mean=2.050947,sd = 1.498367)
        }
        gene_name[i,seq(2,replicate+1)]<-maxt
      }
      for(i in lowBDE){
        maxt<-rnorm(1,mean=2.050947,sd = 1.498367)
        while(maxt>4.536703){
          maxt<-rnorm(1,mean=2.050947,sd = 1.498367)
        }
        gene_name[i,seq(replicate+2,(2*replicate)+1)]<-maxt
      }
      
      temp<-gene_name!=0
      temp[,1]<-F
      gene_name[temp]<-10^as.numeric(gene_name[temp])
      
      for(i in highADE){
        gene_name[i,seq(replicate+2,(2*replicate)+1)]<-gene_name[i,seq(2,replicate+1)]*10^(10^rnorm(n=1,mean=-0.5829726,sd = 0.1976209))
      }
      for(i in highBDE){
        gene_name[i,seq(2,replicate+1)]<-gene_name[i,seq(replicate+2,(2*replicate)+1)]*10^(10^rnorm(n=1,mean=-0.5829726,sd = 0.1976209))
      }
      
      
      for(i in midADE){
        gene_name[i,seq(replicate+2,(2*replicate)+1)]<-gene_name[i,seq(2,replicate+1)]*10^(10^rnorm(n=1,mean=-0.9013155,sd = 0.2863447))
      }
      for(i in midBDE){
        gene_name[i,seq(2,replicate+1)]<-gene_name[i,seq(replicate+2,(2*replicate)+1)]*10^(10^rnorm(n=1,mean=-0.9013155,sd = 0.2863447))
      }
      
      
      for(i in lowADE){
        gene_name[i,seq(replicate+2,(2*replicate)+1)]<-gene_name[i,seq(2,replicate+1)]*10^(10^rnorm(n=1,mean=-1.0660038,sd = 0.4639797))
      }
      for(i in lowBDE){
        gene_name[i,seq(2,replicate+1)]<-gene_name[i,seq(replicate+2,(2*replicate)+1)]*10^(10^rnorm(n=1,mean=-1.0660038,sd = 0.4639797))
      }
      
      for(i in unDE){
        gene_name[i,-1]<-10^rnorm(n = 1,mean = 1.945505,sd = 1.003944)
      }

    
      
    merged_profile<-merge(flux_profile,gene_name,by.x = 'gene_id',by.y = 'gene_name')

##    
    for(row in c(highDE,middleDE,lowDE,unDE)){
      merged_profile[row,group_A]<-merged_profile[row,group_A]+rnorm(length(group_A),mean=0,sd=0.07*sum(merged_profile[row,group_A]))
      merged_profile[row,group_B]<-merged_profile[row,group_B]+rnorm(length(group_B),mean=0,sd=0.07*sum(merged_profile[row,group_B]))
      print(row)
    }
    
    
    for(x in colnames(gene_name[,-1])){
    merged_profile_specified<-data.frame(merged_profile[,1],merged_profile[,2],merged_profile[,3],merged_profile[,4],as.integer(merged_profile[[x]]))
    merged_profile_specified$ratio<-0
    merged_profile_specified$ratio<-sapply(merged_profile_specified[,5], function(y) y/sum(as.numeric(merged_profile_specified[,5]),na.rm = TRUE))
    file<-paste0(name,'/profiled_',x,'.pro')
    write.table(merged_profile_specified[c(2,1,3,4,6,5)],file=file,sep = '\t',row.names=FALSE,quote=FALSE,col.names=F)
    print(x)
    }  
  }
  
  if (!file.exists(name)){
    dir.create(name)
  } 
  
  t<-read.table(profiledata,sep='\t')
  flux_profile<-data.frame(t$V1,t$V2,t$V3,t$V4)
  colnames(flux_profile)<-c('chr_position','gene_id','region','gene_length')
  gene_name<-levels(flux_profile$gene_id)

  set.seed(11111)
  index_gene<-seq(1:length(gene_name))
  highDE<-sample(index_gene,1200)
  middleDE<-sample(index_gene[-highDE],1200)
  DEed<-c(highDE,middleDE)
  lowDE<-sample(index_gene[-DEed],1200)
  DEed<-c(DEed,lowDE)
  unDE<-sample(index_gene[-DEed],2000)


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


  makeprofile()

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