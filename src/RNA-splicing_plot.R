
#summary the numbers of gene occured splicing(all types together) in one project--------
Project_path <- c("ERP107734","SRP011540_Second","SRP070710","SRP094781","SRP150548")

fn_get_count <- function(genes){
  
  data.frame(geneSymbol= unique(genes$geneSymbol),
             Counts=rep(0,length(unique(genes$geneSymbol))))->unique_genes
  
  for (i in 1:nrow(unique_genes)) {
    which(genes$geneSymbol[] == unique_genes$geneSymbol[i]) ->ID
    unique_genes$Counts[i] <- length(ID)
  }
  
  unique_genes
}


for (j in 1:length(Project_path)) {
  list.files(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing",Project_path[j],sep = "/"),
             pattern = "MATS.JCEC.txt",
             recursive = TRUE,full.names = TRUE) ->Project_files
  
  all_genes <- matrix(ncol = 1)%>%as.data.frame(stringsAsFactors=FALSE)
  colnames(all_genes)="geneSymbol"
  
  
  for (i in 1:length(Project_files)) {
    
    read.table(Project_files[i],header = T) %>%
      dplyr::filter(PValue <=0.05)%>%
      #dplyr::filter(IncLevelDifference >= 0.05)%>%
      dplyr::filter(IncLevelDifference <= -0.05)%>%
      dplyr::select(geneSymbol)->sig_symbol
    rbind(sig_symbol,all_genes) -> all_genes
    
  }
  
  fn_get_count(all_genes) -> Project_gene_counts
  Project_gene_counts[order(Project_gene_counts$Counts,decreasing = TRUE),] -> Project_gene_counts
  write.table(Project_gene_counts,col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t",
              paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/",Project_path[j],
                    "/",Project_path[j],"_all_gene_inclusion_counts.txt",sep = ""))#_all_gene_skipping_counts.txt
  
}



# skipping and inclusion seperately-------------------------
rm(list=ls())

#inclusion
list.files("/data/liull/immune-checkpoint-blockade/RNA_splicing",pattern = "_all_gene_inclusion_counts.txt",
           recursive = TRUE,full.names = TRUE)->counts_file
sapply(counts_file, function(x) unlist(strsplit(x,"/"))[6])%>%
  as.vector()->project_id

all_genes_counts=data.frame()
for (i in 1:length(counts_file)) {
  
  read.table(counts_file[i],header = T,sep = "\t",quote = "")%>%
    dplyr::filter(Counts >= 3)->count_file
  count_file$geneSymbol <- as.character(count_file$geneSymbol)
  colnames(count_file)[2] <-  project_id[i]
  
  if(i == 1){
    all_genes_counts <- count_file
  }else{
    all_genes_counts <- merge(all_genes_counts,count_file,all=TRUE)
  }

}

dplyr::mutate(all_genes_counts,SUM=sum(ERP107734 + SRP011540_Second + 
                                         SRP070710 + SRP094781 + SRP150548))->all_genes_counts

dplyr::mutate_all(all_genes_counts,~replace(., is.na(.), 0))%>%
  dplyr::mutate(SUM=ERP107734 + SRP011540_Second + SRP070710 + SRP094781 + SRP150548)->filtered_genes_counts
#%>%dplyr::filter(!SUM == 3)
rownames(filtered_genes_counts) = filtered_genes_counts$geneSymbol
filtered_genes_counts=filtered_genes_counts[,-1]

ID=numeric()
for (i in 1:nrow(filtered_genes_counts)) {
  if(length(which(filtered_genes_counts[i,] == 0)) >= 3){
    ID=c(ID,i)
  }
  
}

filtered_genes_counts[-ID,]->filtered_genes_counts
colnames(filtered_genes_counts)[2]<- "SRP011540"
dplyr::select(filtered_genes_counts,SRP070710,SRP094781,SRP150548,SRP011540,ERP107734)->genes_counts_inclusion


#skipping
list.files("/data/liull/immune-checkpoint-blockade/RNA_splicing",pattern = "_all_gene_skipping_counts.txt",#
           recursive = TRUE,full.names = TRUE)->counts_file
sapply(counts_file, function(x) unlist(strsplit(x,"/"))[6])%>%
  as.vector()->project_id

all_genes_counts=data.frame()
for (i in 1:length(counts_file)) {
  
  read.table(counts_file[i],header = T,sep = "\t",quote = "")%>%
    dplyr::filter(Counts >= 3)->count_file
  count_file$geneSymbol <- as.character(count_file$geneSymbol)
  colnames(count_file)[2] <-  project_id[i]
  
  if(i == 1){
    all_genes_counts <- count_file
  }else{
    all_genes_counts <- merge(all_genes_counts,count_file,all=TRUE)
  }
  
}

dplyr::mutate(all_genes_counts,SUM=sum(ERP107734 + SRP011540_Second + 
                                         SRP070710 + SRP094781 + SRP150548))->all_genes_counts

dplyr::mutate_all(all_genes_counts,~replace(., is.na(.), 0))%>%
  dplyr::mutate(SUM=ERP107734 + SRP011540_Second + SRP070710 + SRP094781 + SRP150548)->filtered_genes_counts
#%>%dplyr::filter(!SUM == 3)
rownames(filtered_genes_counts) = filtered_genes_counts$geneSymbol
filtered_genes_counts=filtered_genes_counts[,-1]

ID=numeric()
for (i in 1:nrow(filtered_genes_counts)) {
  if(length(which(filtered_genes_counts[i,] == 0)) >= 3){
    ID=c(ID,i)
  }
  
}

filtered_genes_counts[-ID,]->filtered_genes_counts
colnames(filtered_genes_counts)[2]<- "SRP011540"
dplyr::select(filtered_genes_counts,SRP070710,SRP094781,SRP150548,SRP011540,ERP107734)->genes_counts_skipping

rm(all_genes_counts,count_file,counts_file,filtered_genes_counts,i,ID,project_id)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


pdf("/data/liull/immune-checkpoint-blockade/RNA_splicing/RNA_splicing_all.pdf",width=6.5,height=5)
my_col_skipping=colorRamp2(c(0,2,4,6),c("#ffffd4","#fed98e","#fe9929","#cc4c02"))#[-c(1,2)]#OrRd
Heatmap(t(genes_counts_skipping),cluster_columns = FALSE,cluster_rows  = FALSE,show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 10),
        col=my_col_skipping,name="Skipping Frequency",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%i", t(genes_counts_skipping)[i, j]), x, y, gp = gpar(fontsize = 10))
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
        })->h1

my_col_inclusion=colorRamp2(c(0,2,4,6,8,10,12),c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594"))#[-c(1,2)]#OrRd
Heatmap(t(genes_counts_inclusion),cluster_columns = FALSE,cluster_rows  = FALSE,
        row_names_gp = gpar(fontsize = 10),
        col=my_col_inclusion,name="Inclusion Frequency",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%i", t(genes_counts_inclusion)[i, j]), x, y, gp = gpar(fontsize = 10))
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
        })->h2
h1 + h2
dev.off()

