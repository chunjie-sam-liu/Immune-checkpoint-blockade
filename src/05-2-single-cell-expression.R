library(tximport)

dir <- "/data/liull/salmon_tutorial/Homo_quants/"
list.files(dir)
sample <- paste0("ERR22287",seq(56,67),"_quant")
sample
files <- file.path(dir,sample,"quant.sf")
all(file.exists(files))

for (i in 1:length(files)) {
  files[i] %>% 
    readr::read_tsv(col_names = TRUE) -> quant
  sapply(quant$Name,function(x) unlist(strsplit(x,"[.]"))[[1]]) %>%
    as.factor() -> quant$Name
  name=strsplit(files[i],"[/]")[[1]][7]
  write.table(quant,paste(dir,name,".sf",sep=""),quote=FALSE,row.names = FALSE,sep="\t")
  
}

files <- file.path(dir,paste(sample,".sf",sep=""))
names(files)=paste0("ERR22287",seq(56,67))

trans_gene = read.table(file.path(dir, "gene_transcript_relationship.txt"),header = T)[,2:1]
txi <- tximport(files, type = "salmon", tx2gene = trans_gene)
write.table(txi$counts,"/data/liull/salmon_tutorial/single_cell_seq_expression.txt",quote=FALSE,sep="\t")

