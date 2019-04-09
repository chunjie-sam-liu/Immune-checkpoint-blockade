library(magrittr)
library(readr)
library(readxl)
library(dplyr)

#filter gastric cancer,RNA-seq,anti-PD1,pretreatment metadata -------------
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="gastric cancer") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

#select response and non_response's sample id and project id-------------------
dplyr::filter(metadata,Response %in% c("CR","PR")) ->response
dplyr::filter(metadata,Response %in% c("SD","PD")) ->non_response

#make difference----------------------------
read.table("/data/liull/immune-checkpoint-blockade/expression/all_count_expression_2.txt",sep="\t",header = T,as.is = TRUE) ->data
all_expression=dplyr::select(data,gene_id,response$Run,non_response$Run)
row.names(all_expression)=all_expression[,1]
all_expression=all_expression[,-1]

DGEList_expr <- DGEList(counts=all_expression)
normalized_expr <- calcNormFactors(DGEList_expr, method="upperquartile")
normalized_loggedCPM_expr = cpm(normalized_expr, log=TRUE, prior.count=2)

keep <- rowSums(normalized_loggedCPM_expr>0) >= 2
normalized_loggedCPM_expr <- normalized_loggedCPM_expr[keep,]
#delete the gene has less than 2 sample exression CPM<1(log2CPM<0)

group_list <- factor(c(rep("response",nrow(response)), rep("non_response",nrow(non_response))))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(normalized_loggedCPM_expr)

fit <- lmFit(normalized_loggedCPM_expr, design)
fit2 <- eBayes(fit)
output <- topTable(fit2, coef=2, n=Inf)
write.table(output,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/all_DEG.txt",quote = FALSE,sep="\t",row.names = TRUE,col.names = TRUE)
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC>1)->up#668
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC< -1)->down#875

######

######

#t.test
# avg.R=apply(normalized_loggedCPM_expr,1,function(x) mean(x[1:nrow(response)]))
# avg.NR=apply(normalized_loggedCPM_expr,1,function(x) mean(x[(nrow(response)+1):ncol(normalized_loggedCPM_expr)]))
# FC=apply(normalized_loggedCPM_expr,1,function(x) (mean(x[1:nrow(response)])+0.01)/(mean(x[(nrow(response)+1):ncol(normalized_loggedCPM_expr)])+0.01))
# p_value=apply(normalized_loggedCPM_expr,1,function(x) t.test(x[1:nrow(response)],x[(nrow(response)+1):ncol(normalized_loggedCPM_expr)])$p.value)
# t_statistic=apply(normalized_loggedCPM_expr,1,function(x) t.test(x[1:nrow(response)],x[(nrow(response)+1):ncol(normalized_loggedCPM_expr)])$statistic)
# result=as.data.frame(cbind(normalized_loggedCPM_expr,avg.R,avg.NR,FC,p_value))
# 
# result=cbind(rownames(result),result)
# rownames(result)=NULL
# colnames(result)[1]="ensembl_ID"
# write.table(result,"/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/gastric_cancer_PD1_DEG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#write all genes' difference
# 
# dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
#   dplyr::filter(FC>=2) -> up#604
# dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
#   dplyr::filter(FC<=0.5) -> down#1102

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,up,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% up$rowname) ->up2
merge(relationship,down,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% down$rowname) ->down2


write.table(up2,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(down2,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#heatmap for all up and down--------------------------------------------------------
rbind(up,down)->all_genes
tibble::rownames_to_column(as.data.frame(normalized_loggedCPM_expr)) %>% 
  dplyr::filter(rowname %in% all_genes$rowname)->expr_heatmap
rownames(expr_heatmap)=expr_heatmap$rowname
expr_heatmap=expr_heatmap[,-1]

apply(expr_heatmap, 1, scale) ->scaled_expr
rownames(scaled_expr)=colnames(expr_heatmap)
scaled_expr=t(scaled_expr)


df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))

pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/heatmap.pdf")
Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 1),col=colorRamp2(c(-4, 0, 4), c("green", "black", "red")))#->origin_heatmap
dev.off()

# > sum(rowSums(scaled_expr>3))
# [1] 140
# > sum(rowSums(scaled_expr< -3))
# [1] 499

# second heatmap

new_scaled_expr <- scaled_expr[row_order(origin_heatmap)[[1]],]

for(i in 1:ncol(new_scaled_expr)) {
  m <- which(new_scaled_expr[,i]>3)
  new_scaled_expr[m,i] <- 3
  n <- which(new_scaled_expr[,i]<(-3))
  new_scaled_expr[n,i] <- (-3)
}
pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/heatmap_2.pdf")
Heatmap(new_scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,cluster_rows = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 1),col=colorRamp2(c(-3, 0, 3), c("green", "black", "red")),heatmap_legend_param=list(at= c(-2, 0, 2)))
dev.off()


#GO enrichment
library(clusterProfiler)
library(org.Hs.eg.db)

enrichGO(gene = up2$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "fdr",pvalueCutoff = 0.05,readable = TRUE)->ego_up#344
DOSE::dotplot(ego_up, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")->ego_up_plot
ggsave(
  filename = 'gastric_cancer_PD1_up_GOenrich.pdf',
  plot = ego_up_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer',
  width = 12,
  height = 8
)

enrichGO(gene = down2$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "fdr",pvalueCutoff = 0.05,readable = TRUE)->ego_down#366
DOSE::dotplot(ego_down, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")->ego_down_plot
ggsave(
  filename = 'gastric_cancer_PD1_down_GOenrich.pdf',
  plot = ego_down_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer',
  width = 12,
  height = 8
)

write.table(as.data.frame(ego_up),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(as.data.frame(ego_down),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)


#KEGG enrichment
enrichKEGG(gene=up2$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH") ->ekegg_up#41
dotplot(ekegg_up, showCategory=20)->KEGG_up_plot
ggsave(
  filename = 'gastric_cancer_PD1_up_KEGG.pdf',
  plot = KEGG_up_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/',
  width = 12,
  height = 8
)
enrichKEGG(gene=down2$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH") ->ekegg_down#21
dotplot(ekegg_down, showCategory=20)->KEGG_down_plot
ggsave(
  filename = 'gastric_cancer_PD1_down_KEGG.pdf',
  plot = KEGG_down_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/',
  width = 12,
  height = 8
)

write.table(as.data.frame(ekegg_up),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(as.data.frame(ekegg_down),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#Reactome enrichment------------------------------------------------------------
library(ReactomePA)
enrichPathway(gene=up2$GeneID,pvalueCutoff=0.05, readable=T)->eReactome_up#33
dotplot(eReactome_up, showCategory=20)->Reactome_up_plot
ggsave(
  filename = 'gastric_cancer_PD1_up_Reactome.pdf',
  plot = Reactome_up_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/',
  width = 12,
  height = 8
)

enrichPathway(gene=down2$GeneID,pvalueCutoff=0.05, readable=T)->eReactome_down#43
dotplot(eReactome_down, showCategory=20)->Reactome_down_plot
ggsave(
  filename = 'gastric_cancer_PD1_down_Reactome.pdf',
  plot = Reactome_down_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/',
  width = 12,
  height = 8
)
write.table(as.data.frame(eReactome_up),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_enrichReactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(as.data.frame(eReactome_down),"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_enrichReactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)


