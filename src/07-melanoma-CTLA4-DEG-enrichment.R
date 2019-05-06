library(magrittr)
library(readr)
library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

#filter melanoma RNA-seq anti-CTLA4 metadata
##use the Second_Response_standard!!!
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Second_Response_standard) ->metadata
metadata %>% dplyr::filter(Run != "SRR3083584") -> metadata# fastq file 16M

dplyr::filter(metadata,Second_Response_standard %in% c("long-survival","R")) -> response
dplyr::filter(metadata,Second_Response_standard %in% c("NR")) -> non_response

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
write.table(output,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/CTLA4_all_DEG.txt",quote = FALSE,row.names = TRUE,col.names = TRUE)

tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<=0.05) %>% dplyr::filter(logFC>1)->up#427
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<=0.05) %>% dplyr::filter(logFC< -1)->down#388
# write.table(up,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_all_up.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)
# write.table(down,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_all_down.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)


read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,up,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% grep("ENSG",up$rowname,value=T)) ->up_ENSG#321
up_ENSG[order(up_ENSG$logFC,decreasing = TRUE),]->up_ENSG
merge(relationship,down,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% grep("ENSG",down$rowname,value=T)) ->down_ENSG#302
down_ENSG[order(down_ENSG$logFC),]->down_ENSG


write.table(up_ENSG,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/CTLA4_up_ENSG.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)
write.table(down_ENSG,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/CTLA4_down_ENSG.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)

#heatmap for ENSG up and down--------------------------------------------------------
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

pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/heatmap_ENSG.pdf")
Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 1),col=colorRamp2(c(-4, 0, 4), c("green", "black", "red")))
dev.off()

# second heatmap
# > sum(rowSums(scaled_expr>2))
# [1] 402
# > sum(rowSums(scaled_expr< -2))
# [1] 774

pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/heatmap_ENSG_2.pdf")
Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 1),col=colorRamp2(c(-2, 0, 2), c("green", "black", "red")))
dev.off()

#GO enrichment
library(clusterProfiler)
library(org.Hs.eg.db)

enrichGO(gene = up_ENSG$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "fdr",pvalueCutoff = 0.05,readable = TRUE)->ego_up#
DOSE::dotplot(ego_up, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")->ego_up_plot
ggsave(
  filename = 'GO_melanoma_CTLA4_up.pdf',
  plot = ego_up_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/',
  width = 12,
  height = 8
)

enrichGO(gene = down_ENSG$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "fdr",pvalueCutoff = 0.05,readable = TRUE)->ego_down#
DOSE::dotplot(ego_down, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")->ego_down_plot
ggsave(
  filename = 'GO_melanoma_CTLA4_down.pdf',
  plot = ego_down_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/',
  width = 12,
  height = 8
)

write.table(as.data.frame(ego_up),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/up_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(as.data.frame(ego_down),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/down_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)


#KEGG enrichment
enrichKEGG(gene=up_ENSG$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH") ->ekegg_up#33
dotplot(ekegg_up, showCategory=20)->KEGG_up_plot
ggsave(
  filename = 'KEGG_melanoma_CTLA4_up.pdf',
  plot = KEGG_up_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/',
  width = 12,
  height = 8
)
enrichKEGG(gene=down_ENSG$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH") ->ekegg_down#5
dotplot(ekegg_down, showCategory=20)->KEGG_down_plot
ggsave(
  filename = 'KEGG_melanoma_CTLA4_down.pdf',
  plot = KEGG_down_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/',
  width = 12,
  height = 8
)

write.table(as.data.frame(ekegg_up),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/up_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(as.data.frame(ekegg_down),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/down_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#Reactome enrichment------------------------------------------------------------
library(ReactomePA)
enrichPathway(gene=up_ENSG$GeneID,pvalueCutoff=0.05, readable=T)->eReactome_up#16
dotplot(eReactome_up, showCategory=20)->Reactome_up_plot
ggsave(
  filename = 'Reactome_melanoma_CTLA4_up.pdf',
  plot = Reactome_up_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/',
  width = 12,
  height = 8
)

enrichPathway(gene=down_ENSG$GeneID,pvalueCutoff=0.05, readable=T)->eReactome_down#4
dotplot(eReactome_down, showCategory=20)->Reactome_down_plot
ggsave(
  filename = 'Reactome_melanoma_CTLA4_down.pdf',
  plot = Reactome_down_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/',
  width = 12,
  height = 8
)

write.table(as.data.frame(eReactome_up),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/up_enrichReactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(as.data.frame(eReactome_down),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/down_enrichReactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)



