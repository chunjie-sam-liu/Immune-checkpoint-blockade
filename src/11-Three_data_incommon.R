# #gene--------------------------------------------------------------------------------------------------------------
# read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_up_ENSG.txt",
#            header = T,as.is = TRUE)->DEG_melanoma_PD1_up
# read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/CTLA4_up_ENSG.txt",
#            header = T,as.is = TRUE)->DEG_melanoma_CTLA4_up
# read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_ENSG.txt",
#            header = T,as.is = TRUE)->DEG_gastric_PD1_up
# venn.diagram(list(melanoma_PD1=melanoma_PD1_up$Ensembl_ID, melanoma_CTLA4=melanoma_CTLA4_up$Ensembl_ID,
#                   gastric_cancer_PD1=gastric_PD1_up$Ensembl_ID),
#              filename=NULL,fill = c("cornflowerblue", "green", "yellow"))->venn_up
# grid.draw(venn_up)
# 
# read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_down_ENSG.txt",
#            header = T,as.is = TRUE)->DEG_melanoma_PD1_down
# read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard//CTLA4_down_ENSG.txt",
#            header = T,as.is = TRUE)->DEG_melanoma_CTLA4_down
# read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_ENSG.txt",
#            header = T,as.is = TRUE)->DEG_gastric_PD1_down
# venn.diagram(list(melanoma_PD1=melanoma_PD1_down$Ensembl_ID, melanoma_CTLA4=melanoma_CTLA4_down$Ensembl_ID,
#                   gastric_cancer_PD1=gastric_PD1_down$Ensembl_ID),
#              filename=NULL,fill = c("cornflowerblue", "green", "yellow"))->venn_down
# grid.draw(venn_down)

#enrich GO-------------------------------------------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/up_enrichGO.txt",
           header = T,as.is = TRUE,sep = "\t",quote = "")%>%
  dplyr::filter(p.adjust <= 0.01)->GO_melanoma_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/up_enrichGO.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->GO_melanoma_CTLA4_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_enrichGO.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->GO_gastric_PD1_up
intersect(GO_melanoma_PD1_up$Description,GO_melanoma_CTLA4_up$Description) %>% intersect(GO_gastric_PD1_up$Description)->up_Intersection

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/down_enrichGO.txt",
           header = T,as.is = TRUE,sep = "\t",quote = "")%>%
  dplyr::filter(p.adjust <= 0.01)->GO_melanoma_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/down_enrichGO.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->GO_melanoma_CTLA4_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_enrichGO.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->GO_gastric_PD1_down
intersect(GO_melanoma_PD1_down$Description,GO_melanoma_CTLA4_down$Description) %>% intersect(GO_gastric_PD1_down$Description)->down_Intersection

#enrichKEGG-----------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/up_enrichKEGG.txt",
           header = T,as.is = TRUE,sep = "\t",quote = "")%>%
  dplyr::filter(p.adjust <= 0.01)->KEGG_melanoma_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/up_enrichKEGG.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->KEGG_melanoma_CTLA4_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_enrichKEGG.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->KEGG_gastric_PD1_up
intersect(KEGG_melanoma_PD1_up$Description,KEGG_melanoma_CTLA4_up$Description) %>% intersect(KEGG_gastric_PD1_up$Description)->KEGG_up_Intersection



read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/down_enrichKEGG.txt",
           header = T,as.is = TRUE,sep = "\t",quote = "")%>%
  dplyr::filter(p.adjust <= 0.01)->KEGG_melanoma_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/down_enrichKEGG.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->KEGG_melanoma_CTLA4_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_enrichKEGG.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->KEGG_gastric_PD1_down
intersect(KEGG_melanoma_PD1_down$Description,KEGG_melanoma_CTLA4_down$Description) %>% intersect(KEGG_gastric_PD1_down$Description)->KEGG_down_Intersection
#0

-----------

fn_all_pathway_genes_GO <- function(enrich_result,GO_Description){
  
  dplyr::filter(enrich_result,Description %in% GO_Description)%>%
    dplyr::select(geneID)-> Genes
  
  sapply(Genes$geneID,function(x) strsplit(x,"/"))%>%
    unlist()%>%
    as.character()%>%
    unique()
  
}
readxl::read_excel("/data/liull/reference/All_EntrezID_Symbl_NCBI.xlsx",col_names = TRUE)%>%as.data.frame(stringsAsFactors)->EntrezID_Symbl
fn_all_pathway_genes_KEGG <- function(enrich_result,KEGG_Description){
  
  dplyr::filter(enrich_result,Description %in% KEGG_Description)%>%
    dplyr::select(geneID)-> Genes
  
  sapply(Genes$geneID,function(x) strsplit(x,"/"))%>%
    unlist()%>%
    as.character()%>%
    unique()->Genes
  data.frame(EntrezID=Genes)%>%
    merge(EntrezID_Symbl,by.x="EntrezID",by.y="GeneID")->Genes
  Genes$Symbol
  
}

#melanoma PD1 map-------
library(ComplexHeatmap)
library(circlize)

m_PD1_up_GO = fn_all_pathway_genes_GO(GO_melanoma_PD1_up,up_Intersection)
m_PD1_down_GO = fn_all_pathway_genes_GO(GO_melanoma_PD1_down,down_Intersection)
m_PD1_up_KEGG = fn_all_pathway_genes_KEGG(KEGG_melanoma_PD1_up,KEGG_up_Intersection)
#m_PD1_down_KEGG:0

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response,Biopsy_Time)%>%
  dplyr::filter(Response != "NE")->metadata
  
dplyr::filter(metadata,Response %in% c("CR","PR","PRCR","R"))-> response
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response


read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_log2CPM_expr.txt",
           header = T)%>%
  dplyr::select(response$Run,non_response$Run)%>%
  tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% c(m_PD1_up_GO,m_PD1_down_GO,m_PD1_up_KEGG))%>%
  tibble::column_to_rownames()->expr_heatmap

# keep <- rowSums(expr_heatmap>0) >= 2
# expr_heatmap <- expr_heatmap[keep,]

apply(expr_heatmap, 1, scale) ->scaled_expr
rownames(scaled_expr)=colnames(expr_heatmap)
scaled_expr=t(scaled_expr)

range(scaled_expr)
sum(rowSums(scaled_expr>2))# 156
sum(rowSums(scaled_expr< -2))#130

df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))

pdf(file="/data/liull/immune-checkpoint-blockade/Three_data_incommon/melanoma_PD1_heatmap.pdf")
Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 3),row_names_gp = gpar(fontsize = 3),col=colorRamp2(c(-2, 0, 2), c("green", "black", "red")))
dev.off()


#melanoma CTLA4 map-------
m_CTLA4_up_GO = fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,up_Intersection)
m_CTLA4_down_GO = fn_all_pathway_genes_GO(GO_melanoma_CTLA4_down,down_Intersection)
m_CTLA4_up_KEGG = fn_all_pathway_genes_KEGG(KEGG_melanoma_CTLA4_up,KEGG_up_Intersection)
#m_CTLA4_down_KEGG 0

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Second_Response_standard) ->metadata
metadata %>% dplyr::filter(Run != "SRR3083584") -> metadata# fastq file 16M

dplyr::filter(metadata,Second_Response_standard %in% c("long-survival","R")) -> response
dplyr::filter(metadata,Second_Response_standard %in% c("NR")) -> non_response

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4_pretreatment_Symbol_count_expr.txt",
           sep="\t",header = T,as.is = TRUE)%>%
  dplyr::select(response$Run,non_response$Run)%>%
  DGEList()%>%
  calcNormFactors(method="upperquartile")%>%
  cpm(log=TRUE, prior.count=2)%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% c(m_CTLA4_up_GO,m_CTLA4_down_GO,m_CTLA4_up_KEGG))%>%
  tibble::column_to_rownames()->expr_heatmap


apply(expr_heatmap, 1, scale) ->scaled_expr
rownames(scaled_expr)=colnames(expr_heatmap)
scaled_expr=t(scaled_expr)

range(scaled_expr)
sum(rowSums(scaled_expr>2))# 
sum(rowSums(scaled_expr< -2))#

df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))

pdf(file="/data/liull/immune-checkpoint-blockade/Three_data_incommon/melanoma_CTLA4_heatmap.pdf")
Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 3),row_names_gp = gpar(fontsize = 3),col=colorRamp2(c(-2, 0, 2), c("green", "black", "red")))
dev.off()


#gastric cancer PD1 map-------
g_PD1_up_GO = fn_all_pathway_genes_GO(GO_gastric_PD1_up,up_Intersection)
g_PD1_down_GO = fn_all_pathway_genes_GO(GO_gastric_PD1_down,down_Intersection)
g_PD1_up_KEGG = fn_all_pathway_genes_KEGG(KEGG_gastric_PD1_up,KEGG_up_Intersection)


readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="gastric cancer") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

dplyr::filter(metadata,Response %in% c("CR","PR")) -> response
dplyr::filter(metadata,Response %in% c("SD","PD")) -> non_response


read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer_PD1_pretreatment_Symbol_count_expr.txt",
           header = T)%>%
  dplyr::select(response$Run,non_response$Run)%>%
  DGEList()%>%
  calcNormFactors(method="upperquartile")%>%
  cpm(log=TRUE, prior.count=2)%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% c(g_PD1_up_GO,g_PD1_down_GO,g_PD1_up_KEGG))%>%
  tibble::column_to_rownames()->expr_heatmap

# keep <- rowSums(expr_heatmap>0) >= 2
# expr_heatmap <- expr_heatmap[keep,]

apply(expr_heatmap, 1, scale) ->scaled_expr
rownames(scaled_expr)=colnames(expr_heatmap)
scaled_expr=t(scaled_expr)

range(scaled_expr)
sum(rowSums(scaled_expr>2))# 
sum(rowSums(scaled_expr< -2))#

df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))

pdf(file="/data/liull/immune-checkpoint-blockade/Three_data_incommon/gastric_cancer_PD1_heatmap.pdf")
Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 3),row_names_gp = gpar(fontsize = 3),col=colorRamp2(c(-2, 0, 2), c("green", "black", "red")))
dev.off()




#ssgaea for gene sets in interaction GO or KEGG pathway-------------------------------------------------
gene_set=list()
gene_set$T_cell_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"T cell activation"),
                                    fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"T cell activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"T cell activation"))
gene_set$regulation_of_T_cell_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of T cell activation"),
                                                  fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of T cell activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of T cell activation"))
gene_set$positive_regulation_of_T_cell_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of T cell activation"),
                                                           fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of T cell activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of T cell activation"))
#
#
gene_set$T_cell_costimulation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"T cell costimulation"),
                                       fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"T cell costimulation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"T cell costimulation"))
#
#
gene_set$lymphocyte_costimulation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"lymphocyte costimulation"),
                                       fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"lymphocyte costimulation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"lymphocyte costimulation"))
gene_set$lymphocyte_differentiation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"lymphocyte differentiation"),
                                           fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"lymphocyte differentiation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"lymphocyte differentiation"))
gene_set$regulation_of_lymphocyte_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of lymphocyte activation"),
                                             fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of lymphocyte activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of lymphocyte activation"))
gene_set$positive_regulation_of_lymphocyte_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of lymphocyte activation"),
                                                               fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of lymphocyte activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of lymphocyte activation"))
gene_set$lymphocyte_proliferation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"lymphocyte proliferation"),
                                           fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"lymphocyte proliferation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"lymphocyte proliferation"))
gene_set$regulation_of_lymphocyte_proliferation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of lymphocyte proliferation"),
                                           fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of lymphocyte proliferation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of lymphocyte proliferation"))
#
#
gene_set$leukocyte_cell_cell_adhesion <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"leukocyte cell-cell adhesion"),
                                           fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"leukocyte cell-cell adhesion"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"leukocyte cell-cell adhesion"))
gene_set$leukocyte_differentiation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"leukocyte differentiation"),
                                           fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"leukocyte differentiation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"leukocyte differentiation"))
gene_set$leukocyte_proliferation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"leukocyte proliferation"),
                                           fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"leukocyte proliferation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"leukocyte proliferation"))
gene_set$positive_regulation_of_leukocyte_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of leukocyte activation"),
                                           fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of leukocyte activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of leukocyte activation"))
gene_set$positive_regulation_of_leukocyte_cell_cell_adhesion <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of leukocyte cell-cell adhesion"),
                                           fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of leukocyte cell-cell adhesion"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of leukocyte cell-cell adhesion"))
gene_set$regulation_of_leukocyte_proliferation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of leukocyte proliferation"),
                                           fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of leukocyte proliferation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of leukocyte proliferation"))
gene_set$regulation_of_leukocyte_cell_cell_adhesion <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of leukocyte cell-cell adhesion"),
                                            fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of leukocyte cell-cell adhesion"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of leukocyte cell-cell adhesion"))
#
#
gene_set$interferon_gamma_production <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"interferon-gamma production"),
                                          fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"interferon-gamma production"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"interferon-gamma production"))
gene_set$regulation_of_interferon_gamma_production <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of interferon-gamma production"),
                                          fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of interferon-gamma production"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of interferon-gamma production"))
#
#
gene_set$mononuclear_cell_proliferation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"mononuclear cell proliferation"),
                                          fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"mononuclear cell proliferation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"mononuclear cell proliferation"))
gene_set$regulation_of_mononuclear_cell_proliferation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of mononuclear cell proliferation"),
                                                     fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of mononuclear cell proliferation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of mononuclear cell proliferation"))
#
#
gene_set$MHC_class_II_protein_complex <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"MHC class II protein complex"),
                                                 fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"MHC class II protein complex"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"MHC class II protein complex"))
#
#
gene_set$regulation_of_cell_cell_adhesion <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of cell-cell adhesion"),
                                             fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of cell-cell adhesion"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of cell-cell adhesion"))
gene_set$positive_regulation_of_cell_adhesion <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of cell adhesion"),
                                                   fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of cell adhesion"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of cell adhesion"))
gene_set$positive_regulation_of_cell_cell_adhesion <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of cell-cell adhesion"),
                                             fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of cell-cell adhesion"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of cell-cell adhesion"))
#
gene_set$antigen_receptor_mediated_signaling_pathway <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"antigen receptor-mediated signaling pathway"),
                                                 fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"antigen receptor-mediated signaling pathway"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"antigen receptor-mediated signaling pathway"))
gene_set$cellular_defense_response <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"cellular defense response"),
                                                 fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"cellular defense response"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"cellular defense response"))
gene_set$immune_response_activating_cell_surface_receptor_signaling_pathway <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"immune response-activating cell surface receptor signaling pathway"),
                                                 fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"immune response-activating cell surface receptor signaling pathway"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"immune response-activating cell surface receptor signaling pathway"))
gene_set$immune_response_regulating_cell_surface_receptor_signaling_pathway <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"immune response-regulating cell surface receptor signaling pathway"),
                                                 fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"immune response-regulating cell surface receptor signaling pathway"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"immune response-regulating cell surface receptor signaling pathway"))
gene_set$immunological_synapse <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"immunological synapse"),
                                                 fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"immunological synapse"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"immunological synapse"))
gene_set$positive_regulation_of_cell_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of cell activation"),
                                        fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of cell activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of cell activation"))
gene_set$regulation_of_immune_effector_process <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of immune effector process"),
                                        fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of immune effector process"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of immune effector process"))
gene_set$response_to_lipopolysaccharide <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"response to lipopolysaccharide"),
                                        fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"response to lipopolysaccharide"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"response to lipopolysaccharide"))
gene_set$response_to_molecule_of_bacterial_origin <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"response to molecule of bacterial origin"),
                                        fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"response to molecule of bacterial origin"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"response to molecule of bacterial origin"))
gene_set$side_of_membrane <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"side of membrane"),
                                        fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"side of membrane"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"side of membrane"))

#
#
gene_set$Hematopoietic_cell_lineage <- union(fn_all_pathway_genes_KEGG(KEGG_melanoma_PD1_up,"Hematopoietic cell lineage"),
                                             fn_all_pathway_genes_KEGG(KEGG_melanoma_CTLA4_up,"Hematopoietic cell lineage"))%>%
  union(fn_all_pathway_genes_KEGG(KEGG_gastric_PD1_up,"Hematopoietic cell lineage"))
gene_set$Staphylococcus_aureus_infection <- union(fn_all_pathway_genes_KEGG(KEGG_melanoma_PD1_up,"Staphylococcus aureus infection"),
                                             fn_all_pathway_genes_KEGG(KEGG_melanoma_CTLA4_up,"Staphylococcus aureus infection"))%>%
  union(fn_all_pathway_genes_KEGG(KEGG_gastric_PD1_up,"Staphylococcus aureus infection"))
gene_set$Intestinal_immune_network_for_IgA_production <- union(fn_all_pathway_genes_KEGG(KEGG_melanoma_PD1_up,"Intestinal immune network for IgA production"),
                                             fn_all_pathway_genes_KEGG(KEGG_melanoma_CTLA4_up,"Intestinal immune network for IgA production"))%>%
  union(fn_all_pathway_genes_KEGG(KEGG_gastric_PD1_up,"Intestinal immune network for IgA production"))
gene_set$CAMs <- union(fn_all_pathway_genes_KEGG(KEGG_melanoma_PD1_up,"Cell adhesion molecules (CAMs)"),
                                             fn_all_pathway_genes_KEGG(KEGG_melanoma_CTLA4_up,"Cell adhesion molecules (CAMs)"))%>%
  union(fn_all_pathway_genes_KEGG(KEGG_gastric_PD1_up,"Cell adhesion molecules (CAMs)"))
gene_set$Th17_cell_differentiation <- union(fn_all_pathway_genes_KEGG(KEGG_melanoma_PD1_up,"Th17 cell differentiation"),
                                             fn_all_pathway_genes_KEGG(KEGG_melanoma_CTLA4_up,"Th17 cell differentiation"))%>%
  union(fn_all_pathway_genes_KEGG(KEGG_gastric_PD1_up,"Th17 cell differentiation"))


#melanoma PD1 ssgsea test
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_log2CPM_expr.txt",
           header = T)->melanoma_PD1_log2CPM

ssgava_score <- gsva(as.matrix(melanoma_PD1_log2CPM), gene_set, min.sz=1, max.sz=999999, method="ssgsea",kcdf="Gaussian")


readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "SRA")%>%
  dplyr::filter(Cancer == "melanoma")%>%
  dplyr::filter(Anti_target =="anti-PD1")%>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(Run,SRA_Study,Response)%>%
  dplyr::filter(Response != "NE")->metadata
metadata$Response %>%
  gsub("PD","NR",.)%>% gsub("SD","NR",.)%>%
  gsub("PR","R",.)%>%gsub("CR","R",.)->metadata$Response
dplyr::filter(metadata,Response=="R") ->R
dplyr::filter(metadata,Response=="NR")->NR
as.data.frame(ssgava_score,stringsAsFactors=FALSE)%>%
  dplyr::select(R$Run,NR$Run)->ordered_ssgava_score

for (i in 1:nrow(ordered_ssgava_score)) {
  
  as.numeric(ordered_ssgava_score[i,1:nrow(R)])->R_numer
  as.numeric(ordered_ssgava_score[i,(nrow(R)+1):ncol(ordered_ssgava_score)])->NR_numer
  wilcox.test(R_numer,NR_numer)$p.value->P
  if(P< 0.05){
    print(P)
    print(rownames(ordered_ssgava_score)[i])
  }
  
  
}


#melanoma CTLA4 ssgsea test
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4_pretreatment_Symbol_count_expr.txt",
           header = T)->melanoma_CTLA4_count
DGEList_expr <- DGEList(counts=melanoma_CTLA4_count)
normalized_expr <- calcNormFactors(DGEList_expr, method="upperquartile")
melanoma_CTLA4_log2CPM = cpm(normalized_expr, log=TRUE, prior.count=2)
write.table(melanoma_CTLA4_log2CPM,
            "/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4_pretreatment_Symbol_log2CPM_expr.txt",
            quote = FALSE,row.names = TRUE,col.names = TRUE)


ssgava_score <- gsva(as.matrix(melanoma_CTLA4_log2CPM), gene_set, min.sz=1, max.sz=999999, method="ssgsea",kcdf="Gaussian")


readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "dbGAP")%>%
  dplyr::filter(Cancer == "melanoma")%>%
  dplyr::filter(Anti_target =="anti-CTLA4")%>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(Run,SRA_Study,Second_Response_standard)%>%
  dplyr::filter(Run != "SRR3083584")->metadata
metadata$Second_Response_standard %>%
  gsub("long-survival","R",.)->metadata$Second_Response_standard
dplyr::filter(metadata,Second_Response_standard=="R") ->R
dplyr::filter(metadata,Second_Response_standard=="NR")->NR
as.data.frame(ssgava_score,stringsAsFactors=FALSE)%>%
  dplyr::select(R$Run,NR$Run)->ordered_ssgava_score

for (i in 1:nrow(ordered_ssgava_score)) {
  
  as.numeric(ordered_ssgava_score[i,1:nrow(R)])->R_numer
  as.numeric(ordered_ssgava_score[i,(nrow(R)+1):ncol(ordered_ssgava_score)])->NR_numer
  wilcox.test(R_numer,NR_numer)$p.value->P
  if(P< 0.05){
    print(P)
    print(rownames(ordered_ssgava_score)[i])
  }
  
  
}


#gastric PD1 ssgsea test
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer_PD1_pretreatment_Symbol_count_expr.txt",
           header = T)->gastric_PD1_count
DGEList_expr <- DGEList(counts=gastric_PD1_count)
normalized_expr <- calcNormFactors(DGEList_expr, method="upperquartile")
gastric_PD1_log2CPM = cpm(normalized_expr, log=TRUE, prior.count=2)
write.table(gastric_PD1_log2CPM,
            "/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer_PD1_pretreatment_Symbol_log2CPM_expr.txt",
            quote = FALSE,row.names = TRUE,col.names = TRUE)

ssgava_score <- gsva(as.matrix(gastric_PD1_log2CPM), gene_set, min.sz=1, max.sz=999999, method="ssgsea",kcdf="Gaussian")


readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "SRA")%>%
  dplyr::filter(Cancer == "gastric cancer")%>%
  dplyr::filter(Anti_target =="anti-PD1")%>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(Run,SRA_Study,Response)->metadata
metadata$Response %>%
  gsub("PD","NR",.)%>% gsub("SD","NR",.)%>%
  gsub("PR","R",.)%>%gsub("CR","R",.)->metadata$Response
dplyr::filter(metadata,Response=="R") ->R
dplyr::filter(metadata,Response=="NR")->NR
as.data.frame(ssgava_score,stringsAsFactors=FALSE)%>%
  dplyr::select(R$Run,NR$Run)->ordered_ssgava_score

for (i in 1:nrow(ordered_ssgava_score)) {
  
  as.numeric(ordered_ssgava_score[i,1:nrow(R)])->R_numer
  as.numeric(ordered_ssgava_score[i,(nrow(R)+1):ncol(ordered_ssgava_score)])->NR_numer
  wilcox.test(R_numer,NR_numer)$p.value->P
  if(P< 0.05){
    print(P)
    print(rownames(ordered_ssgava_score)[i])
  }
  
  
}




