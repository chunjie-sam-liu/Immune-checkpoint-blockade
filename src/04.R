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


assem_count_files %>%
  head(1) %>%
  read_tsv(col_names = TRUE) %>%
  as.data.frame() -> reads_file
  ID=grep("_mRNA",reads_file$accession)    
  mRNA_reads_file=reads_file[ID,]
      
      
      
      
      
       
 