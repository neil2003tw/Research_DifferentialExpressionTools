
recursivecomparison<-function(normalnum,controlnum,data_folder,output_folder='./'){


#source('/NFSShare/neil/DEtools/working_dir/DifferentialExpressionTools/R scripts/Comparison/Tools/Odds.R')
#source('/NFSShare/neil/DEtools/working_dir/DifferentialExpressionTools/R scripts/Comparison/Tools/Pipe.R')
#source('/NFSShare/neil/DEtools/working_dir//DifferentialExpressionTools/R scripts/Comparison/Tools/DEcomparison.R')
  
  
source('/Users/NeilWu/Github/DifferentialExpressionTools/R scripts/Comparison/Tools/Odds.R')
source('/Users/NeilWu/Github/DifferentialExpressionTools/R scripts/Comparison/Tools/Pipe.R')
source('/Users/NeilWu/Github/DifferentialExpressionTools/R scripts/Comparison/Tools/DEcomparison.R')


count_table<-counttable_merged(data_folder)
AnswerDE<-read.csv(paste0(data_folder,'/allDElist.csv'))

### adjust Answer matrix
AnswerDE<-AnswerDE[!is.na(AnswerDE$DE),]
levels(AnswerDE$DE)<-c('highDE','highDE','lowDE','lowDE','midDE','midDE','unDE')

### move to working folder
origindir<-getwd()
setwd(output_folder)

mergeddata<-data.frame(DEcomparison(count_table,AnswerDE,normalnum,controlnum))
mergeddata$Replicate_count<-normalnum
mergeddata$DEtools<-row.names(mergeddata)

for(x in seq(normalnum-1,1)){
  nextdata<-data.frame(DEcomparison(count_table,AnswerDE,x,x))
  nextdata$Replicate_count<-x
  nextdata$DEtools<-row.names(nextdata)
  mergeddata<-rbind(mergeddata,nextdata)
}

### back to original folder
setwd(origindir)

return(mergeddata)
}
