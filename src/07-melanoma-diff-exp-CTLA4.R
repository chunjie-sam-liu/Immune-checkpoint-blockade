library(magrittr)
library(readr)
library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

#filter melanoma RNA-seq anti-CTLA4 metadata
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata
metadata %>% dplyr::filter(Run != "SRR3083584") -> metadata# fastq file 16M

dplyr::filter(metadata,Response %in% c("CR","PR","X")) -> response
dplyr::filter(metadata,Response %in% c("SD","PD")) -> non_response

#load the origin counts data
read.table("/data/liull/immune-checkpoint-blockade/expression/all_count_expression_2.txt",sep="\t",header = T,as.is = TRUE) ->all_counts_data
dplyr::select(all_counts_data,gene_id,response$Run,non_response$Run)->ordered_melanoma_CTLA4
row.names(ordered_melanoma_CTLA4)=ordered_melanoma_CTLA4[,1]
ordered_melanoma_CTLA4=ordered_melanoma_CTLA4[,-1]

DGEList_expr <- DGEList(counts=ordered_melanoma_CTLA4)
normalized_expr <- calcNormFactors(DGEList_expr, method="upperquartile")
normalized_loggedCPM_expr = cpm(normalized_expr, log=TRUE, prior.count=2)
dim(normalized_loggedCPM_expr)

keep <- rowSums(normalized_loggedCPM_expr>0) >= 2
normalized_loggedCPM_expr <- normalized_loggedCPM_expr[keep,]
dim(normalized_loggedCPM_expr)
#delete the gene has less than 2 sample exression CPM<1(log2CPM<0)

group_list <- factor(c(rep("response",nrow(response)), rep("non_response",nrow(non_response))))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(normalized_loggedCPM_expr)

fit <- lmFit(normalized_loggedCPM_expr, design)
fit2 <- eBayes(fit)
output <- topTable(fit2, coef=2, n=Inf)

tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<=0.05) %>% dplyr::filter(logFC>log2(1.5))->up#273
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<=0.05) %>% dplyr::filter(logFC< -log2(1.5))->down#270

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI_2.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,up,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% grep("ENSG",up$rowname,value=T)) ->up_ENSG#200
up_ENSG[order(up_ENSG$logFC,decreasing = TRUE),]->up_ENSG
merge(relationship,down,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% grep("ENSG",down$rowname,value=T)) ->down_ENSG#187
down_ENSG[order(down_ENSG$logFC),]->down_ENSG

# write.table(output,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/all_DEG.txt",quote = FALSE,sep="\t",row.names = TRUE,col.names = TRUE)
# write.table(up_ENSG,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_ENSG.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)
# write.table(down_ENSG,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_ENSG.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)

#heatmap for all up and down--------------------------------------------------------
rbind(up_ENSG,down_ENSG)->all_genes
tibble::rownames_to_column(as.data.frame(normalized_loggedCPM_expr)) %>% 
  dplyr::filter(rowname %in% all_genes$Ensembl_ID)->expr_heatmap
rownames(expr_heatmap)=expr_heatmap$rowname
expr_heatmap=expr_heatmap[,-1]

apply(expr_heatmap, 1, scale) ->scaled_expr
rownames(scaled_expr)=colnames(expr_heatmap)
scaled_expr=t(scaled_expr)


df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))

pdf(file="/data/liull/immune-checkpoint-blockade/heatmap_ENSG.pdf")
Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 1),col=colorRamp2(c(-4, 0, 4), c("green", "black", "red")))
dev.off()

# second heatmap
# > sum(rowSums(scaled_expr>2))
# [1] 293
# > sum(rowSums(scaled_expr< -2))
# [1] 353

pdf(file="/data/liull/immune-checkpoint-blockade/heatmap_ENSG_2.pdf")
Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 1),col=colorRamp2(c(-2, 0, 2), c("green", "black", "red")))
dev.off()

#GO enrichment
library(clusterProfiler)
library(org.Hs.eg.db)

enrichGO(gene = up_ENSG$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "fdr",pvalueCutoff = 0.05,readable = TRUE)->ego_up#7
DOSE::dotplot(ego_up, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")->ego_up_plot
ggsave(
  filename = 'melanoma_CTLA4_up_GOenrich.pdf',
  plot = ego_up_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/',
  width = 12,
  height = 8
)

enrichGO(gene = down_ENSG$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "fdr",pvalueCutoff = 0.05,readable = TRUE)->ego_down#13
DOSE::dotplot(ego_down, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")->ego_down_plot
ggsave(
  filename = 'melanoma_CTLA4_down_GOenrich.pdf',
  plot = ego_down_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/',
  width = 12,
  height = 8
)

# write.table(as.data.frame(ego_up),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
# write.table(as.data.frame(ego_down),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)


#KEGG enrichment
enrichKEGG(gene=up_ENSG$GeneID,organism="human",pvalueCutoff=1,pAdjustMethod = "BH") ->ekegg_up#0
# dotplot(ekegg_up, showCategory=20)->KEGG_up_plot
# ggsave(
#   filename = 'gastric_cancer_PD1_up_KEGG.pdf',
#   plot = KEGG_up_plot,
#   device = 'pdf',
#   path = '/data/liull/immune-checkpoint-blockade/',
#   width = 12,
#   height = 8
# )
enrichKEGG(gene=down_ENSG$GeneID,organism="human",pvalueCutoff=1,pAdjustMethod = "BH") ->ekegg_down#0
# dotplot(ekegg_down, showCategory=20)->KEGG_down_plot
# ggsave(
#   filename = 'gastric_cancer_PD1_down_KEGG.pdf',
#   plot = KEGG_down_plot,
#   device = 'pdf',
#   path = '/data/liull/immune-checkpoint-blockade/',
#   width = 12,
#   height = 8
# )

# write.table(as.data.frame(ekegg_up),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
# write.table(as.data.frame(ekegg_down),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#Reactome enrichment------------------------------------------------------------
library(ReactomePA)
enrichPathway(gene=up_ENSG$GeneID,pvalueCutoff=1, readable=T)->eReactome_up#HDMs demethylate histones; p.adjust=0.07799137
# dotplot(eReactome_up, showCategory=20)->Reactome_up_plot
# ggsave(
#   filename = 'gastric_cancer_PD1_up_Reactome.pdf',
#   plot = Reactome_up_plot,
#   device = 'pdf',
#   path = '/data/liull/immune-checkpoint-blockade/',
#   width = 12,
#   height = 8
# )

enrichPathway(gene=down_ENSG$GeneID,pvalueCutoff=1, readable=T)->eReactome_down#Cell surface interactions at the vascular wall; p.adjust=0.02255345 
# dotplot(eReactome_down, showCategory=20)->Reactome_down_plot
# ggsave(
#   filename = 'gastric_cancer_PD1_down_Reactome.pdf',
#   plot = Reactome_down_plot,
#   device = 'pdf',
#   path = '/data/liull/immune-checkpoint-blockade/',
#   width = 12,
#   height = 8
# )
# write.table(as.data.frame(eReactome_up),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_enrichReactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
# write.table(as.data.frame(eReactome_down),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_enrichReactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)



