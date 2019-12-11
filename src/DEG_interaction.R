#DEG gene venn plot---------------------------------------------------------------------
library (VennDiagram)

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_up_ENSG.txt",header = T)%>%
  dplyr::select(Symbol)%>%
  dplyr::filter(Symbol != "NA")->melanoma_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_down_ENSG.txt",header = T)%>%
  dplyr::select(Symbol)%>%
  dplyr::filter(Symbol != "NA")->melanoma_PD1_down

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_ENSG.txt",header = T)%>%
  dplyr::select(Symbol)%>%
  dplyr::filter(Symbol != "NA")->gastric_cancer_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_ENSG.txt",header = T)%>%
  dplyr::select(Symbol)%>%
  dplyr::filter(Symbol != "NA")->gastric_cancer_down

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/CTLA4_up_ENSG.txt",
           header = T)%>%dplyr::select(Symbol)%>%
  dplyr::filter(Symbol != "NA")->melanoma_CTLA4_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/CTLA4_down_ENSG.txt",
           header = T)%>%dplyr::select(Symbol)%>%
  dplyr::filter(Symbol != "NA")->melanoma_CTLA4_down

# intersect(melanoma_PD1_up$Symbol,gastric_cancer_up$Symbol) %>% intersect(melanoma_CTLA4_up$Symbol)->all_up_genes
# #"UBD"      "LCK"      "JAKMIP1"  "TRAT1"    "IDO1"     "HLA-DOA"  "PDCD1"    "HLA-DQA1"
# intersect(melanoma_PD1_down$Symbol,gastric_cancer_down$Symbol) %>% intersect(melanoma_CTLA4_down$Symbol)->all_down_genes
# #"CORO2B"

venn.diagram(list(gastric_PD1=gastric_cancer_up$Symbol, 
                  melanoma_PD1=melanoma_PD1_up$Symbol,
                  melanoma_CTLA4=melanoma_CTLA4_up$Symbol),
             filename=NULL,fill = c("#66c2a5", "#fc8d62", "#8da0cb"),
             alpha = 0.90,cex = 2.0,cat.cex = 1.5,cat.dist = 0.1)->venn_up
#alpha:颜色透明度；cex：数字大小；cat.cex：类别名字大小；cat.dist：类别名距离图像的距离
pdf("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/venn_up_genes.pdf")
grid.draw(venn_up)
dev.off()



venn.diagram(list(gastric_PD1=gastric_cancer_down$Symbol, 
                  melanoma_PD1=melanoma_PD1_down$Symbol,
                  melanoma_CTLA4=melanoma_CTLA4_down$Symbol),
             filename=NULL,fill = c("#66c2a5", "#fc8d62", "#8da0cb"),
             alpha = 0.90,cex = 2.0,cat.cex = 1.5,cat.dist = 0.1)->venn_down
pdf("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/venn_down_genes.pdf")
grid.draw(venn_down)
dev.off()
#


# #intersection genes' heatmap plot (single)-----------------------------------------------------------------------------------------------
# all_genes <- c("UBD","LCK","JAKMIP1","TRAT1","IDO1","HLA-DOA","PDCD1","HLA-DQA1","CORO2B")
# 
# 
# #melanoma PD1
# readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
#   dplyr::filter(Library_strategy=="RNA-Seq") %>%
#   dplyr::filter(Cancer=="melanoma") %>%
#   dplyr::filter(Anti_target=="anti-PD1") %>%
#   dplyr::filter(Biopsy_Time == "pre-treatment")%>%
#   dplyr::select(Run,Response) ->metadata
# dplyr::filter(metadata,Response %in% c("CR","PR","R"))-> response
# dplyr::filter(metadata,Response %in% c("SD","PD","NR"))-> non_response
# 
# read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_log2CPM_expr.txt",
#            sep = "\t",header = T)%>%tibble::rownames_to_column()%>%
#   dplyr::filter(rowname %in% all_genes)%>%
#   tibble::column_to_rownames()%>%
#   dplyr::select(response$Run,non_response$Run)->melanoma_PD1_log2CPM_expr
# apply(melanoma_PD1_log2CPM_expr, 1, scale) ->scaled_expr
# rownames(scaled_expr)=colnames(melanoma_PD1_log2CPM_expr)
# scaled_expr=t(scaled_expr)
# 
# sum(rowSums(scaled_expr>2))
# sum(rowSums(scaled_expr< -2))
# 
# df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
# ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))
# #mycol=colorRamp2(c(2,1,0,-1,-2),brewer.pal(11,"RdYlBu")[-c(1,11)])
# mycol=colorRamp2(c(2,0,-2),c("#E56663","white","#4488C9"))
# 
# pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_9_heatmap.pdf",width = 18,height = 8)
# Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,cluster_rows = FALSE,
#         column_names_gp = gpar(fontsize = 10),row_names_gp = gpar(fontsize = 15),
#         col=mycol)
# dev.off()
# 
# #melanoma CTLA4
# readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP") %>%
#   dplyr::filter(Library_strategy=="RNA-Seq") %>%
#   dplyr::filter(Cancer=="melanoma") %>%
#   dplyr::filter(Anti_target=="anti-CTLA4") %>%
#   dplyr::filter(Biopsy_Time == "pre-treatment")%>%
#   dplyr::select(Run,Second_Response_standard) ->metadata
# metadata %>% dplyr::filter(Run != "SRR3083584") -> metadata# fastq file 16M
# dplyr::filter(metadata,Second_Response_standard %in% c("long-survival","R"))-> response
# dplyr::filter(metadata,Second_Response_standard %in% c("NR"))-> non_response
# 
# read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4_pretreatment_Symbol_log2CPM_expr.txt",
#            header = T)%>%tibble::rownames_to_column()%>%
#   dplyr::filter(rowname %in% all_genes)%>%
#   tibble::column_to_rownames()%>%
#   dplyr::select(response$Run,non_response$Run)->melanoma_CTLA4_log2CPM_expr
# apply(melanoma_CTLA4_log2CPM_expr, 1, scale) ->scaled_expr
# rownames(scaled_expr)=colnames(melanoma_CTLA4_log2CPM_expr)
# scaled_expr=t(scaled_expr)
# 
# range(scaled_expr)
# sum(rowSums(scaled_expr>2))
# sum(rowSums(scaled_expr< -2))
# 
# df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
# ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))
# #mycol=colorRamp2(c(2,1,0,-1,-2),brewer.pal(11,"RdYlBu")[-c(1,11)])
# mycol=colorRamp2(c(2,0,-2),c("#E56663","white","#4488C9"))
# 
# pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4_9_heatmap.pdf",width = 18,height = 8)
# Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,cluster_rows = FALSE,
#         column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 15),
#         col=mycol)
# dev.off()
# 
# #gastric PD1
# readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
#   dplyr::filter(Library_strategy=="RNA-Seq") %>%
#   dplyr::filter(Cancer=="gastric cancer") %>%
#   dplyr::filter(Anti_target=="anti-PD1") %>%
#   dplyr::filter(Biopsy_Time == "pre-treatment")%>%
#   dplyr::select(Run,Response) ->metadata
# dplyr::filter(metadata,Response %in% c("CR","PR"))-> response
# dplyr::filter(metadata,Response %in% c("SD","PD"))-> non_response
# 
# read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer_PD1_pretreatment_Symbol_log2CPM_expr.txt",
#            header = T)%>%tibble::rownames_to_column()%>%
#   dplyr::filter(rowname %in% all_genes)%>%
#   tibble::column_to_rownames()%>%
#   dplyr::select(response$Run,non_response$Run)->gastric_PD1_log2CPM_expr
# apply(gastric_PD1_log2CPM_expr, 1, scale) ->scaled_expr
# rownames(scaled_expr)=colnames(gastric_PD1_log2CPM_expr)
# scaled_expr=t(scaled_expr)
# 
# range(scaled_expr)
# sum(rowSums(scaled_expr>2))
# sum(rowSums(scaled_expr< -2))
# 
# df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
# ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))
# #mycol=colorRamp2(c(2,1,0,-1,-2),brewer.pal(11,"RdYlBu")[-c(1,11)])
# mycol=colorRamp2(c(2,0,-2),c("#E56663","white","#4488C9"))
# 
# pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer_PD1_9_heatmap.pdf",width = 18,height = 8)
# Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,cluster_rows = FALSE,
#         column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 15),
#         col=mycol)
# dev.off()


#intersection genes' heatmap plot (all together)-----------------------------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
all_genes <- c("UBD","LCK","JAKMIP1","TRAT1","IDO1","HLA-DOA","PDCD1","HLA-DQA1","CORO2B")

#melanoma PD1
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(Run,Response) ->metadata
dplyr::filter(metadata,Response %in% c("CR","PR","R"))-> response_1
dplyr::filter(metadata,Response %in% c("SD","PD","NR"))-> non_response_1

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_log2CPM_expr.txt",
           sep = "\t",header = T)%>%tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% all_genes)%>%
  tibble::column_to_rownames()->melanoma_PD1
colnames(melanoma_PD1) -> melanoma_PD1_samples
apply(melanoma_PD1,1, scale)%>%t()%>%as.data.frame() ->melanoma_PD1
colnames(melanoma_PD1) <- melanoma_PD1_samples

dplyr::select(melanoma_PD1,response_1$Run)%>%
  apply(., 1,mean)->melanoma_aPD1_R
dplyr::select(melanoma_PD1,non_response_1$Run)%>%
  apply(., 1,mean)->melanoma_aPD1_NR

cbind(melanoma_aPD1_R,melanoma_aPD1_NR)->expr

#melanoma CTLA4
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(Run,Second_Response_standard) ->metadata
metadata %>% dplyr::filter(Run != "SRR3083584") -> metadata# fastq file 16M
dplyr::filter(metadata,Second_Response_standard %in% c("long-survival","R"))-> response_2
dplyr::filter(metadata,Second_Response_standard %in% c("NR"))-> non_response_2

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4_pretreatment_Symbol_log2CPM_expr.txt",
           header = T)%>%tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% all_genes)%>%
  tibble::column_to_rownames()->melanoma_CTLA4
colnames(melanoma_CTLA4) -> melanoma_CTLA4_samples
apply(melanoma_CTLA4,1, scale)%>%t()%>%as.data.frame() ->melanoma_CTLA4
colnames(melanoma_CTLA4) <- melanoma_CTLA4_samples

dplyr::select(melanoma_CTLA4,response_2$Run)%>%
  apply(., 1,mean)->melanoma_aCTLA4_R
dplyr::select(melanoma_CTLA4,non_response_2$Run)%>%
  apply(., 1,mean)->melanoma_aCTLA4_NR

cbind(expr,melanoma_aCTLA4_R)%>%cbind(melanoma_aCTLA4_NR)->expr


#gastric
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="gastric cancer") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(Run,Response) ->metadata
dplyr::filter(metadata,Response %in% c("CR","PR"))-> response_3
dplyr::filter(metadata,Response %in% c("SD","PD"))-> non_response_3

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer_PD1_pretreatment_Symbol_log2CPM_expr.txt",
           header = T)%>%tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% all_genes)%>%
  tibble::column_to_rownames()->gastric_PD1
colnames(gastric_PD1) -> gastric_PD1_samples
apply(gastric_PD1,1, scale)%>%t()%>%as.data.frame() ->gastric_PD1
colnames(gastric_PD1) <- gastric_PD1_samples

dplyr::select(gastric_PD1,response_3$Run)%>%
  apply(., 1,mean)->gastric_aPD1_R
dplyr::select(gastric_PD1,non_response_3$Run)%>%
  apply(., 1,mean)->gastric_aPD1_NR

cbind(expr,gastric_aPD1_R)%>%cbind(gastric_aPD1_NR)->expr

#


range(expr)
# sum(rowSums(expr>2))
# sum(rowSums(expr< -2))

# df = data.frame(type = c(rep("response", nrow(response_1)), rep("non_response", nrow(non_response_1)),
#                          rep("response", nrow(response_2)), rep("non_response", nrow(non_response_2)),
#                          rep("response", nrow(response_3)), rep("non_response", nrow(non_response_3))))
# ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))
#mycol=colorRamp2(c(2,1,0,-1,-2),brewer.pal(11,"RdYlBu")[-c(1,11)])
mycol=colorRamp2(c(0.8,0,-0.8),c("#E56663","white","#4488C9"))
# data.frame(Class=c(rep("up",8),"down"))->Split
# rownames(Split) <- all_genes

pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/All_9_heatmap.pdf")
Heatmap(expr,name="Color_key",cluster_columns = FALSE,cluster_rows = FALSE,
        column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 10),
        col=mycol,heatmap_legend_param = list(border="black",color_bar="continous"))
dev.off()
