#melanoma PD1 removed batch effect expr--------------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/PD1_removed_batch_expression.txt",
           header = T)->melanoma_PD1
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE)->relationship

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response,Biopsy_Time) ->metadata

tibble::rownames_to_column(melanoma_PD1)%>%
  dplyr::select(rowname,metadata$Run)%>%
  dplyr::filter(rowname %in% grep("ENSG",rownames(melanoma_PD1),value=T)) ->ENSG_melanoma_PD1

relationship%>%  
  dplyr::filter(Ensembl_ID %in% ENSG_melanoma_PD1$rowname)%>%
  dplyr::select(Symbol,Ensembl_ID)%>%
  merge(ENSG_melanoma_PD1,by.x="Ensembl_ID",by.y = "rowname")%>%
  dplyr::select(-Ensembl_ID)->Symbol_melanoma_PD1

factors=factor(Symbol_melanoma_PD1$Symbol)
merged_expression=tapply(Symbol_melanoma_PD1[,2],factors,median)
for (i in 3:ncol(Symbol_melanoma_PD1)) {
  temp=tapply(Symbol_melanoma_PD1[,i],factors,median)
  merged_expression=cbind(merged_expression,temp)
}
colnames(merged_expression)=colnames(Symbol_melanoma_PD1)[2:ncol(Symbol_melanoma_PD1)]
dim(merged_expression)
write.table(merged_expression,
            "/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")

#melanoma CTLA4 expr--------------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/expression/all_FPKM_expression_2.txt",
           header = T)->all_FPKM_expr
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE)->relationship

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Second_Response_standard) ->metadata
metadata %>% dplyr::filter(Run != "SRR3083584") -> metadata# fastq file 16M

all_FPKM_expr%>%
  dplyr::select(gene_id,metadata$Run)%>%
  dplyr::filter(gene_id %in% grep("ENSG",all_FPKM_expr$gene_id,value=T)) ->ENSG_melanoma_CTLA4

relationship%>%  
  dplyr::filter(Ensembl_ID %in% ENSG_melanoma_CTLA4$gene_id)%>%
  dplyr::select(Symbol,Ensembl_ID)%>%
  merge(ENSG_melanoma_CTLA4,by.x="Ensembl_ID",by.y = "gene_id")%>%
  dplyr::select(-Ensembl_ID)->Symbol_melanoma_CTLA4

factors=factor(Symbol_melanoma_CTLA4$Symbol)
merged_expression=tapply(Symbol_melanoma_CTLA4[,2],factors,median)
for (i in 3:ncol(Symbol_melanoma_CTLA4)) {
  temp=tapply(Symbol_melanoma_CTLA4[,i],factors,median)
  merged_expression=cbind(merged_expression,temp)
}
colnames(merged_expression)=colnames(Symbol_melanoma_CTLA4)[2:ncol(Symbol_melanoma_CTLA4)]
dim(merged_expression)
write.table(merged_expression,
            "/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4_pretreatment_Symbol_FPKM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")


#gastric cancer PD1 expr--------------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/expression/all_FPKM_expression_2.txt",
           header = T)->all_FPKM_expr
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE)->relationship

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="gastric cancer") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Second_Response_standard) ->metadata

all_FPKM_expr%>%
  dplyr::select(gene_id,metadata$Run)%>%
  dplyr::filter(gene_id %in% grep("ENSG",all_FPKM_expr$gene_id,value=T)) ->ENSG_gastric_cancer_PD1

relationship%>%  
  dplyr::filter(Ensembl_ID %in% ENSG_gastric_cancer_PD1$gene_id)%>%
  dplyr::select(Symbol,Ensembl_ID)%>%
  merge(ENSG_gastric_cancer_PD1,by.x="Ensembl_ID",by.y = "gene_id")%>%
  dplyr::select(-Ensembl_ID)->Symbol_gastric_cancer_PD1

factors=factor(Symbol_gastric_cancer_PD1$Symbol)
merged_expression=tapply(Symbol_gastric_cancer_PD1[,2],factors,median)
for (i in 3:ncol(Symbol_gastric_cancer_PD1)) {
  temp=tapply(Symbol_gastric_cancer_PD1[,i],factors,median)
  merged_expression=cbind(merged_expression,temp)
}
colnames(merged_expression)=colnames(Symbol_gastric_cancer_PD1)[2:ncol(Symbol_gastric_cancer_PD1)]
dim(merged_expression)
write.table(merged_expression,
            "/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer_PD1_pretreatment_Symbol_FPKM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")

#melanoma PD1 two project(not remove batch effect)---------------------
read.table("/data/liull/immune-checkpoint-blockade/expression/all_FPKM_expression_2.txt",
           header = T)->all_FPKM_expr
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE)->relationship

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Second_Response_standard) ->metadata

all_FPKM_expr%>%
  dplyr::select(gene_id,metadata$Run)%>%
  dplyr::filter(gene_id %in% grep("ENSG",all_FPKM_expr$gene_id,value=T)) ->ENSG_melanoma_PD1

relationship%>%  
  dplyr::filter(Ensembl_ID %in% ENSG_melanoma_PD1$gene_id)%>%
  dplyr::select(Symbol,Ensembl_ID)%>%
  merge(ENSG_melanoma_PD1,by.x="Ensembl_ID",by.y = "gene_id")%>%
  dplyr::select(-Ensembl_ID)->Symbol_melanoma_PD1

factors=factor(Symbol_melanoma_PD1$Symbol)
merged_expression=tapply(Symbol_melanoma_PD1[,2],factors,median)
for (i in 3:ncol(Symbol_melanoma_PD1)) {
  temp=tapply(Symbol_melanoma_PD1[,i],factors,median)
  merged_expression=cbind(merged_expression,temp)
}
colnames(merged_expression)=colnames(Symbol_melanoma_PD1)[2:ncol(Symbol_melanoma_PD1)]
dim(merged_expression)
write.table(merged_expression,
            "/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_FPKM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")
