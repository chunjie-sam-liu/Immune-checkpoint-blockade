library(magrittr)
library(readr)
library(readxl)
library(dplyr)
library(clusterProfiler)
library(ReactomePA)

#load melanoma_PD1_removed_batch profile
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_PD1_removed_batch_expression.txt",header = T,as.is = TRUE)->combat_edata

#filter melanoma RNA-seq anti-PD1's response information
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%###
  dplyr::select(SRA_Study,Run,Response) ->metadata

#select and order profile--response~nonresponse
dplyr::filter(metadata,Response %in% c("CR","PR","PRCR","R")) -> response#26
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response#59
dplyr::select(combat_edata,response$Run,non_response$Run)->ordered_melanoma_PD1

avg.R=apply(ordered_melanoma_PD1,1,function(x) median(x[1:nrow(response)]))
avg.NR=apply(ordered_melanoma_PD1,1,function(x) median(x[(nrow(response)+1):length(ordered_melanoma_PD1)]))
diff.avg=apply(ordered_melanoma_PD1,1,function(x) (median(x[1:nrow(response)])-median(x[(nrow(response)+1):length(ordered_melanoma_PD1)])))
p_value=apply(ordered_melanoma_PD1,1,function(x) t.test(x[1:nrow(response)],x[(nrow(response)+1):length(ordered_melanoma_PD1)])$p.value)
t_statistic=apply(ordered_melanoma_PD1,1,function(x) t.test(x[1:nrow(response)],x[(nrow(response)+1):length(ordered_melanoma_PD1)])$statistic)

result=as.data.frame(cbind(ordered_melanoma_PD1,avg.R,avg.NR,diff.avg,p_value,t_statistic))
result=cbind(rownames(result),result)
rownames(result)=NULL
colnames(result)[1]="ensembl_ID"
write.table(result,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/pretreatment/melanoma_PD1_DEG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#write all genes' difference


dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
  dplyr::filter(diff.avg>=2) -> up#36

dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
  dplyr::filter(diff.avg<=-2) -> down#78

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,up,by.x="Ensembl_ID",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% up$ensembl_ID) ->up2
merge(relationship,down,by.x="Ensembl_ID",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% down$ensembl_ID) ->down2
write.table(up2,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/pretreatment/up.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(down2,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/pretreatment/down.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)


#GO enrichment--------------------------------------------------------------------------------------
enrichGO(gene = up2$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,readable = TRUE)->ego_up
dotplot(ego_up,showCategory=20)->ego_up_plot
enrichGO(gene = down2$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,readable = TRUE) ->ego_down
dotplot(ego_down,showCategory=20)->ego_down_plot

ggsave(
  filename = 'melanoma_PD1_down_GOenrich.pdf',
  plot = ego_down_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/pretreatment/',
  width = 12,
  height = 8
)

write.table(as.data.frame(ego_up),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/pretreatment/up_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#
write.table(as.data.frame(ego_down),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/pretreatment/down_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#

#KEGG enrichment
enrichKEGG(gene=up2$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH") ->ekegg_up#32
enrichKEGG(gene=down2$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH")->ekegg_down#1
browseKEGG(ekegg_up, 'hsa04934')
browseKEGG(ekegg_down, 'hsa04934')
write.table(as.data.frame(ekegg_up),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/up_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(as.data.frame(ekegg_down),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/down_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#Reactome enrichment------------------------------------------------------
enrichPathway(gene=up2$GeneID,organism="human",pvalueCutoff=0.05, readable=T)->eReactome_up
dotplot(eReactome_up, showCategory=20)#0
enrichPathway(gene=down2$GeneID,organism="human",pvalueCutoff=0.05, readable=T)->eReactome_down
dotplot(eReactome_down, showCategory=20)->Reactome_down_plot#13

ggsave(
  filename = 'melanoma_PD1_down_Reactome.pdf',
  plot = Reactome_down_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/pretreatment/',
  width = 12,
  height = 8
)

write.table(as.data.frame(eReactome_up),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/pretreatment/up_enrichReactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(as.data.frame(eReactome_down),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/pretreatment/down_enrichReactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)




#Heatmap(ComplexHeatmap)
rbind(up,down)->all_genes
rownames(all_genes)=all_genes$ensembl_ID
all_genes[,-c(1,(ncol(all_genes)-4),(ncol(all_genes)-3),(ncol(all_genes)-2),(ncol(all_genes)-1),(ncol(all_genes)))]->all_genes

apply(all_genes, 1, scale) ->scared_all_genes
rownames(scared_all_genes)=colnames(all_genes)
scared_all_genes=t(scared_all_genes)

Heatmap(scared_all_genes,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 3),row_names_gp = gpar(fontsize = 3))


#add labels
rbind(as.matrix(c(rep("response",nrow(response)),rep("non_response",nrow(non_response))))[,1],all_genes) %>%
  t()%>%
  as.data.frame()->labels
annot_up <- data.frame(response_VS_none = labels$`1`)
colors = list(response_VS_none = c("response" = "red", "non_response" = "blue"))
ha <- HeatmapAnnotation(annot_up, col = colors)

pdf(file="/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/pretreatment/heatmap.pdf")
Heatmap(scared_all_genes,column_names_gp = gpar(fontsize = 3),row_names_gp = gpar(fontsize = 3),top_annotation = ha,cluster_columns = FALSE)
dev.off()

# #pheatmap
# 
# rbind(up,down)->all_genes
# rownames(all_genes)=all_genes$ensembl_ID
# all_genes[,-c(1,(ncol(all_genes)-4),(ncol(all_genes)-3),(ncol(all_genes)-2),(ncol(all_genes)-1),(ncol(all_genes)))]->all_genes
# annotation_col = data.frame(GeneClass = factor(rep(c("response", "non-response"), c(nrow(response),nrow(non_response)))))
# rownames(annotation_col)=colnames(all_genes)
# pheatmap(all_genes, annotation_col = annotation_col,cluster_cols = FALSE,scale="row")

