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
            "/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_log2CPM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_log2CPM_expr.txt",
           header = T)->melanoma_PD1_pretreatment_Symbol_log2CPM_expr

metadata%>%
  dplyr::filter(SRA_Study == "SRP070710")->SRP070710_metadata
SRP070710_metadata %>%
  dplyr::filter(Response %in% c("CR","PR"))->SRP070710_R
SRP070710_metadata %>%
  dplyr::filter(Response %in% c("SD","PD"))->SRP070710_NR
melanoma_PD1_pretreatment_Symbol_log2CPM_expr %>%
  dplyr::select(SRP070710_metadata$Run)->SRP070710_log2CPM_expr
melanoma_PD1_pretreatment_Symbol_log2CPM_expr %>%
  dplyr::select(SRP070710_R$Run)->SRP070710_Symbol_log2CPM_expr_R
melanoma_PD1_pretreatment_Symbol_log2CPM_expr %>%
  dplyr::select(SRP070710_NR$Run)->SRP070710_Symbol_log2CPM_expr_NR

write.table(SRP070710_log2CPM_expr,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/SRP070710_Symbol_log2CPM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")
write.table(SRP070710_Symbol_log2CPM_expr_R,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP070710_log2CPM/SRP070710_Symbol_logCPM_expr_R.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")
write.table(SRP070710_Symbol_log2CPM_expr_NR,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_NR/SRP070710_log2CPM/SRP070710_Symbol_logCPM_expr_NR.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")


metadata%>%
  dplyr::filter(SRA_Study == "SRP094781")->SRP094781_metadata
SRP094781_metadata %>%
  dplyr::filter(Response %in% c("CR","PR"))->SRP094781_R
SRP094781_metadata %>%
  dplyr::filter(Response %in% c("SD","PD"))->SRP094781_NR

melanoma_PD1_pretreatment_Symbol_log2CPM_expr %>%
  dplyr::select(SRP094781_metadata$Run)->SRP094781_log2CPM_expr
melanoma_PD1_pretreatment_Symbol_log2CPM_expr %>%
  dplyr::select(SRP094781_R$Run)->SRP094781_Symbol_log2CPM_expr_R
melanoma_PD1_pretreatment_Symbol_log2CPM_expr %>%
  dplyr::select(SRP094781_NR$Run)->SRP094781_Symbol_log2CPM_expr_NR

write.table(SRP094781_log2CPM_expr,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/SRP094781_log2CPM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")
write.table(SRP094781_Symbol_log2CPM_expr_R,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781_log2CPM/SRP094781_Symbol_log2CPM_expr_R.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")
write.table(SRP094781_Symbol_log2CPM_expr_NR,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_NR/SRP094781_log2CPM/SRP094781_Symbol_log2CPM_expr_NR.txt",
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

#melanoma PD1 two project FPKM(not remove batch effect)---------------------
read.table("/data/liull/immune-checkpoint-blockade/expression/all_FPKM_expression_2.txt",
           header = T)->all_FPKM_expr
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE)->relationship

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

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

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_FPKM_expr.txt",
           header = T)->melanoma_PD1_pretreatment_Symbol_FPKM_expr
metadata%>%
  dplyr::filter(SRA_Study == "SRP070710")->SRP070710_metadata
SRP070710_metadata %>%
  dplyr::filter(Response %in% c("CR","PR"))->SRP070710_R
SRP070710_metadata %>%
  dplyr::filter(Response %in% c("SD","PD"))->SRP070710_NR
melanoma_PD1_pretreatment_Symbol_FPKM_expr %>%
  dplyr::select(SRP070710_metadata$Run)->SRP070710_expr
melanoma_PD1_pretreatment_Symbol_FPKM_expr %>%
  dplyr::select(SRP070710_R$Run)->SRP070710_Symbol_FPKM_expr_R
melanoma_PD1_pretreatment_Symbol_FPKM_expr %>%
  dplyr::select(SRP070710_NR$Run)->SRP070710_Symbol_FPKM_expr_NR

write.table(SRP070710_expr,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/SRP070710_Symbol_FPKM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")
write.table(SRP070710_Symbol_FPKM_expr_R,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP070710/SRP070710_Symbol_FPKM_expr_R.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")
write.table(SRP070710_Symbol_FPKM_expr_NR,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_NR/SRP070710/SRP070710_Symbol_FPKM_expr_NR.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")

metadata%>%
  dplyr::filter(SRA_Study == "SRP094781")->SRP094781_metadata
SRP094781_metadata %>%
  dplyr::filter(Response %in% c("CR","PR"))->SRP094781_R
SRP094781_metadata %>%
  dplyr::filter(Response %in% c("SD","PD"))->SRP094781_NR

melanoma_PD1_pretreatment_Symbol_FPKM_expr %>%
  dplyr::select(SRP094781_metadata$Run)->SRP094781_expr
melanoma_PD1_pretreatment_Symbol_FPKM_expr %>%
  dplyr::select(SRP094781_R$Run)->SRP094781_Symbol_FPKM_expr_R
melanoma_PD1_pretreatment_Symbol_FPKM_expr %>%
  dplyr::select(SRP094781_NR$Run)->SRP094781_Symbol_FPKM_expr_NR

write.table(SRP094781_expr,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/SRP094781_Symbol_FPKM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")
write.table(SRP094781_Symbol_FPKM_expr_R,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/SRP094781_Symbol_FPKM_expr_R.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")
write.table(SRP094781_Symbol_FPKM_expr_NR,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_NR/SRP094781/SRP094781_Symbol_FPKM_expr_NR.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")

#melanoma PD1 two project count(not remove batch effect)---------------------
read.table("/data/liull/immune-checkpoint-blockade/expression/all_count_expression_2.txt",
           header = T)->all_count_expr
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE)->relationship

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Second_Response_standard) ->metadata

all_count_expr%>%
  dplyr::select(gene_id,metadata$Run)%>%
  dplyr::filter(gene_id %in% grep("ENSG",all_count_expr$gene_id,value=T)) ->ENSG_melanoma_PD1

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
            "/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_count_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")


#single project of FPKM,pretreatment,symbol-----------------------------------
read.table("/data/liull/immune-checkpoint-blockade/expression/all_FPKM_expression_2.txt",
           header = T)->all_FPKM_expr
read.table("/data/liull/reference/RNA_seq_pipeline_EnsemblID_Symbl.txt",header = T,as.is = TRUE)->relationship


#melanoma CTLA4:SRP011540
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Second_Response_standard)%>% 
  dplyr::filter(Run != "SRR3083584")-> metadata# fastq file 16M
#melanoma PD1:SRP070710
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(SRA_Study == "SRP070710")%>%
  dplyr::select(SRA_Study,Run) ->metadata
#melanoma PD1:SRP094781
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(SRA_Study == "SRP094781")%>%
  dplyr::select(SRA_Study,Run) ->metadata
#melanoma PD1:SRP150548
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(SRA_Study == "SRP150548")%>%
  dplyr::select(SRA_Study,Run) ->metadata

#gastric cancer PD1:ERP107734
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="gastric cancer") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(SRA_Study == "ERP107734")%>%
  dplyr::select(SRA_Study,Run) ->metadata


all_FPKM_expr%>%
  dplyr::select("gene_id",metadata$Run)%>%
  dplyr::filter(gene_id %in% grep("ENSG",all_FPKM_expr$gene_id,value=T))->ENSG_expr

relationship%>%  
  dplyr::filter(Ensembl_ID %in% ENSG_expr$gene_id)%>%
  dplyr::select(Symbol,Ensembl_ID)%>%
  merge(ENSG_expr,by.x="Ensembl_ID",by.y = "gene_id")%>%
  dplyr::select(-Ensembl_ID)->Symbol_expr

factors=factor(Symbol_expr$Symbol)
merged_expr=tapply(Symbol_expr[,2],factors,median)
for (i in 3:ncol(Symbol_expr)) {
  temp=tapply(Symbol_expr[,i],factors,median)
  merged_expr=cbind(merged_expr,temp)
}
colnames(merged_expr)=colnames(Symbol_expr)[2:ncol(Symbol_expr)]
dim(merged_expr)

write.table(merged_expr,
            "/data/liull/immune-checkpoint-blockade/expression/SRP011540_pretreatment_Symbol_FPKM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")

write.table(merged_expr,
            "/data/liull/immune-checkpoint-blockade/expression/SRP070710_pretreatment_Symbol_FPKM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")
write.table(merged_expr,
            "/data/liull/immune-checkpoint-blockade/expression/SRP094781_pretreatment_Symbol_FPKM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")
write.table(merged_expr,
            "/data/liull/immune-checkpoint-blockade/expression/SRP150548_pretreatment_Symbol_FPKM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")

write.table(merged_expr,
            "/data/liull/immune-checkpoint-blockade/expression/ERP107734_pretreatment_Symbol_FPKM_expr.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")




