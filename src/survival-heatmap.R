library(magrittr)
library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)


read.table("/data/liull/immune-checkpoint-blockade/survival/PD1_survival_symbol.txt",header = T,as.is = TRUE)->PD1_symbol
read.table("/data/liull/immune-checkpoint-blockade/survival/CTLA4_survival_symbol.txt",header = T,as.is = TRUE)->CTLA4_symbol
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_PD1_removed_batch_expression.txt",header = T,as.is = TRUE)->PD1_expr
cbind(rownames(PD1_expr),PD1_expr)->PD1_expr
colnames(PD1_expr)[1]="Ensembl_ID"
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_CTLA4_removed_batch_expression.txt",header = T,as.is = TRUE)->CTLA4_expr
cbind(rownames(CTLA4_expr),CTLA4_expr)->CTLA4_expr
colnames(CTLA4_expr)[1]="Ensembl_ID"

read.table("/data/liull/reference/EnsemblID_Symbl_ensembl.txt",sep="\t",header = T,as.is = TRUE) ->relationship

#filter melanoma RNA-seq response information
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") ->SRA
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP") ->dbGAP

#PD1------------------------------------------------------------------------------------------------
SRA %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%###
  dplyr::select(SRA_Study,Run,Response) ->PD1_metadata

#select and order profile--response~nonresponse
dplyr::filter(PD1_metadata,Response %in% c("CR","PR","PRCR","R")) -> PD1_response#26
dplyr::filter(PD1_metadata,Response %in% c("SD","PD","NR")) -> PD1_non_response#59

#
relationship %>%
  merge(PD1_expr)%>%
  dplyr::filter(Symbol %in% PD1_symbol$gene_symbol)%>%
  dplyr::select(Symbol,PD1_response$Run,PD1_non_response$Run)->PD1_symbol_expr

rownames(PD1_symbol_expr)=PD1_symbol_expr$Symbol
PD1_symbol_expr[,-1]->PD1_symbol_expr
annotation_col = data.frame(SampleClass = factor(rep(c("response", "non-response"), c(nrow(PD1_response),nrow(PD1_non_response)))))
rownames(annotation_col)=colnames(PD1_symbol_expr)
pheatmap(PD1_symbol_expr, annotation_col = annotation_col,scale="row")->PD1_survival_plot
ggsave(
  filename = 'PD1_survival_heatmap.pdf',
  plot = PD1_survival_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/survival/',
  width = 16,
  height = 6.8
)


#CTLA4---------------------------------------------------------------
SRA %>%
  rbind(dbGAP)%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%###
  dplyr::select(SRA_Study,Run,Response) ->CTLA4_metadata

#select and order profile--response~nonresponse
dplyr::filter(CTLA4_metadata,Response %in% c("CR","PR","PRCR","R")) -> CTLA4_response#37
dplyr::filter(CTLA4_metadata,Response %in% c("SD","PD","NR")) -> CTLA4_non_response#7

#
relationship %>%
  merge(CTLA4_expr)%>%
  dplyr::filter(Symbol %in% CTLA4_symbol$gene_symbol)%>%
  dplyr::select(Symbol,CTLA4_response$Run,CTLA4_non_response$Run)->CTLA4_symbol_expr

rownames(CTLA4_symbol_expr)=CTLA4_symbol_expr$Symbol
CTLA4_symbol_expr[,-1]->CTLA4_symbol_expr
annotation_col = data.frame(SampleClass = factor(rep(c("response", "non-response"), c(nrow(CTLA4_response),nrow(CTLA4_non_response)))))
rownames(annotation_col)=colnames(CTLA4_symbol_expr)
pheatmap(CTLA4_symbol_expr, annotation_col = annotation_col,scale="row")->CTLA4_survival_plot#cluster_cols = FALSE
ggsave(
  filename = 'CTLA4_survival_heatmap.pdf',
  plot = CTLA4_survival_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/survival/',
  width = 16,
  height = 6.8
)
