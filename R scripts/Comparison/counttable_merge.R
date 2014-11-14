
##DElist<-DElist[!is.na(DElist$DE),]
##test[is.na(test)]<-0
counttable_merged<-function(directory){
  
setwd(directory)

data_name<-dir()
data_name<-data_name[grepl('profiled',data_name)]
condition<-substr(data_name,10,11)
datas<-paste0('profile',condition)
for(x in seq(1,length(condition))){
  assign(datas[x],read.csv(data_name[x],sep = ' ',col.names=c('Gene',condition[x])))
}
##list_profile<-lapply(ls(), function(x) if (class(get(x)) == "data.frame") get(x))
list_profile<-lapply(ls(pattern = 'profile[AB][0-9]'), function(x) get(x))
count_table<-Reduce(function(x, y) merge(x, y, all=TRUE), list_profile)
count_table[is.na(count_table)]<-0
count_table<-data.frame(count_table[,-1],row.names=count_table[,1])
rm(list=datas)
setwd('../R scripts/Comparison/')

return(count_table)
}
