#melanoma PD1 removed batch effect expr
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/PD1_removed_batch_expression.txt",
           header = T)->melanoma_PD1
read.table("/data/liull/reference/EnsemblID_Symbl_ensembl.txt",header = T,as.is = TRUE)->relationship

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
  merge(ENSG_melanoma_PD1,by.x="Ensembl_ID",by.y = "rowname",all = TRUE)%>%
  dplyr::select(-Ensembl_ID)->Symbol_melanoma_PD1

factors=factor(Symbol_melanoma_PD1$Symbol)
merged_expression=tapply(Symbol_melanoma_PD1[,2],factors,median)
for (i in 3:ncol(Symbol_melanoma_PD1)) {
  temp=tapply(Symbol_melanoma_PD1[,i],factors,median)
  merged_expression=cbind(merged_expression,temp)
}
colnames(merged_expression)=colnames(Symbol_melanoma_PD1)[2:ncol(Symbol_melanoma_PD1)]




