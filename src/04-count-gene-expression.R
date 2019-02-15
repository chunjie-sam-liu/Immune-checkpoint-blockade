# library -----------------------------------------------------------------

library(magrittr)
library(dplyr)
library(purrr)
# path --------------------------------------------------------------------

path_data <- '/data/liucj/data/immune-checkpoint-blockade/fastq/'

# load gtf expression -----------------------------------------------------

assem_count_files <- list.files(path = path_data, pattern = 'read_count.xls', recursive = TRUE, full.names = TRUE)[-c(1:12)]

# get gene and transcript id relationship from ----------------------------

readr::read_tsv("/project/huhui/lnc_pipeline_v2/LNC_pipeline_DB/Homo/Anno/all_rmdup_RNA_Ensembl_Chr_Id_GRCh38_81_seleno_filted.gtf",col_names = FALSE, skip = 5) %>% 
  as.data.frame() %>%
  dplyr::filter(X3 =='transcript') -> .gtf
mrna = grep("ENST",.gtf$X9,value=T)
.gtf %>%
  dplyr::filter(X9 %in% mrna) ->mrna_relation
as.vector(sapply(mrna_relation$X9,function(x) unlist(strsplit(x,"\""))[2])) -> gene_ID
as.vector(sapply(mrna_relation$X9,function(x) unlist(strsplit(x,"\""))[6])) -> transcript_ID
relationship1=cbind(gene_ID,transcript_ID)
ncrna = grep("NONH",.gtf$X9,value=T)
.gtf %>%
  dplyr::filter(X9 %in% ncrna) ->ncrna_relation
as.vector(sapply(ncrna_relation$X9,function(x) unlist(strsplit(x,"\""))[2])) -> gene_ID
as.vector(sapply(ncrna_relation$X9,function(x) unlist(strsplit(x,"\""))[4])) -> transcript_ID
relationship2=cbind(gene_ID,transcript_ID)
relationship=rbind(relationship1,relationship2)
write.table(relationship,"/data/liull/immune-checkpoint-blockade/RNA_seq_pipeline_gtf_relationship.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)


all_expression=matrix(NA)
colnames(all_expression)="gene_id"
#
for (i in 1:length(assem_count_files)) {
assem_count_files[i] %>%
  read_tsv(col_names = TRUE) %>%
  as.data.frame() -> reads_file

#filter the transcript id ended with "_mRNA" and delete the "_mRNA"
ID=grep("_mRNA",reads_file$accession)    
mRNA_reads_file=reads_file[ID,]
as.vector(sapply(mRNA_reads_file$accession,function(x) unlist(strsplit(x,"_"))[1])) ->transcript_ID    
mRNA_reads_file$accession = transcript_ID
#compute the all count number
dplyr::mutate(mRNA_reads_file,all_count=tag_count_Forward+tag_count_Reverse) -> mRNA_reads_file

#add the ensembl gene id 
expression=merge(mRNA_reads_file,relationship,by.x="accession",by.y="transcript_ID")

##the same one gene's count value added ,result in a new file only with gene id and count number
expression2=as.data.frame(tapply(expression$all_count,factor(expression$gene_ID),sum))      
expression2=cbind(row.names(expression2), expression2)
row.names(expression2)=NULL 
name=strsplit(strsplit(assem_count_files[i],"/")[[1]][14],"_")[[1]][1]
colnames(expression2)=c("gene_id",name)      

#merge with others
all_expression=merge(all_expression,expression2,by="gene_id",all=T)

}
write.table(all_expression,"/data/liull/immune-checkpoint-blockade/expression/all_count_expression_2.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

 