
##DElist<-DElist[!is.na(DElist$DE),]
##test[is.na(test)]<-0

data_name<-dir()
condition<-substr(dir(),10,11)
datas<-paste0('profile',condition)
for(name in seq(1,length(dir()))){
  assign(datas[name],read.csv(data_name[name],sep = ' ',col.names=c('Gene',condition[name])))
}
##list_profile<-lapply(ls(), function(x) if (class(get(x)) == "data.frame") get(x))
list_profile<-lapply(ls(pattern = 'profile[AB][0-9]'), function(x) get(x))
count_table<-Reduce(function(x, y) merge(x, y, all=TRUE), list_profile)
count_table[is.na(count_table)]<-0
count_table<-data.frame(count_table[,-1],row.names=count_table[,1])
