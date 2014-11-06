count_fluxdata_fastq<-function(data){
 ###Reading fastq and extract transcript id
 awk_data_name<-paste0("awk 'BEGIN{i=0}{i++;if (i%4==1) print $1}' < ",data)
 file.pipe <- pipe(awk_data_name)
 fastq_line1<-read.table(file.pipe)
 sep_fastq_line1<-strsplit(as.character(fastq_line1$V1),':')
 transcript_id_raw<-sapply(sep_fastq_line1,function(x) x[3])
 
 ###Creating transcript_id to gene_name matrix
 #annotation<-read.table(annotation_file,sep='\t')
 #gene_name_ann<-sapply(strsplit(as.character(annotation$V9),';'),'[[',2)
 #gene_name_ann<-sapply(strsplit(gene_name_ann,' '),'[[',3)
 #transcript_id_ann<-sapply(strsplit(as.character(annotation$V9),';'),'[[',1)
 #transcript_id_ann<-sapply(strsplit(transcript_id_ann,' '),'[[',2)
 #annotation<-data.frame(transcript_id_ann,gene_name_ann)
 #annotation<-annotation[!duplicated(annotation),]
 
 ###Convert
 table()
 
 data_name<-paste0(data,"_counted.txt")
 write.table(table(factor(transcript_id_raw)),file=data_name,row.names=FALSE,quote=FALSE)
}

dir_list<-dir()
dir_list<-dir_list[grep('.*./fastq',dir_list)]
apply(dir_list,2,count_fluxdata_fastq)
