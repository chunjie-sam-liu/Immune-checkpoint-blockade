
# library -----------------------------------------------------------------

library(magrittr)
library(dplyr)

# path --------------------------------------------------------------------

path_data <- '/data/liucj/data/immune-checkpoint-blockade/fastq/'

# load gtf expression -----------------------------------------------------

assem_gtf_files <- list.files(path = path_data, pattern = 'assem.gtf', recursive = TRUE, full.names = TRUE)
#select transcript----------------------------------------------------------------------------------

all_rna_expression=matrix(NA)
colnames(all_rna_expression)="gene_id"


for (i in 13:length(assem_gtf_files)) {
  

assem_gtf_files[i] %>% 
  readr::read_tsv(col_names = FALSE, skip = 2)  %>% 
  dplyr::filter(X3 =='transcript') -> transcript

#select coding gene----------------------------------------------------------------------------------
mrna = grep("ENST",transcript$X9,value=T)
transcript %>%
  dplyr::filter(X9 %in% mrna) ->mrna_expression

as.vector(sapply(mrna_expression$X9,function(x) unlist(strsplit(x,";"))[[3]])) %>% 
  sapply(function(x) unlist(strsplit(x,"\""))[[2]]) %>%
  as.vector() -> reference_id

as.vector(sapply(mrna_expression$X9,function(x) unlist(strsplit(x,";"))[[4]])) %>% 
  sapply(function(x) unlist(strsplit(x,"\""))[[2]]) %>%
  as.vector() ->ref_gene_id

as.vector(sapply(mrna_expression$X9,function(x) unlist(strsplit(x,";"))[[5]])) %>% 
  sapply(function(x) unlist(strsplit(x,"\""))[[2]]) %>%
  as.vector() ->ref_gene_name

as.vector(sapply(mrna_expression$X9,function(x) unlist(strsplit(x,";"))[[7]])) %>% 
  sapply(function(x) unlist(strsplit(x,"\""))[[2]]) %>%
  as.vector() ->FPKM
cbind(reference_id,ref_gene_id,ref_gene_name,FPKM) -> mrna_FPKM
#select non-coding gene-----------------------------------------------------------------------------
ncrna = grep("NONHSAT",transcript$X9,value=T)
transcript %>%
  dplyr::filter(X9 %in% ncrna) ->ncrna_expression

as.vector(sapply(ncrna_expression$X9,function(x) unlist(strsplit(x,";"))[[3]])) %>% 
  sapply(function(x) unlist(strsplit(x,"\""))[[2]]) %>%
  as.vector() -> ncrna_reference_id

as.vector(sapply(ncrna_expression$X9,function(x) unlist(strsplit(x,";"))[[4]])) %>% 
  sapply(function(x) unlist(strsplit(x,"\""))[[2]]) %>%
  as.vector() ->ncrna_ref_gene_id

ncrna_ref_gene_name = rep(NA,length(ncrna_ref_gene_id))

as.vector(sapply(ncrna_expression$X9,function(x) unlist(strsplit(x,";"))[[6]])) %>% 
  sapply(function(x) unlist(strsplit(x,"\""))[[2]]) %>%
  as.vector() ->nc_FPKM_value
cbind(ncrna_reference_id,ncrna_ref_gene_id,ncrna_ref_gene_name,nc_FPKM_value) -> ncrna_FPKM

all=rbind(mrna_FPKM,ncrna_FPKM)

#the same one gene's FPKM value added-----------------------------------------------------------------
all_expression=as.data.frame(tapply(as.numeric(all[,4]),factor(all[,2]),sum))

all_expression=cbind(row.names(all_expression), all_expression)
row.names(all_expression)=NULL
name=strsplit(strsplit(assem_gtf_files[i],"/")[[1]][14],"\\.")[[1]][1]
colnames(all_expression)=c("gene_id",name)

#所有样本表达值组合
all_rna_expression=merge(all_rna_expression,all_expression,by="gene_id",all=T)

}

write.table(all_rna_expression[1:(dim(all_rna_expression)[1]-1),],"/data/liull/immune-checkpoint-blockade/expression/all_FPKM_expression_2.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

