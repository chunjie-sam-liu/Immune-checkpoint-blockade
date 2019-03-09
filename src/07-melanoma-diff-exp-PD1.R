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
  dplyr::filter(diff.avg>=2) -> up

dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
  dplyr::filter(diff.avg<=-2) -> down

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,up,by.x="EnsemblId",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(EnsemblId %in% up$ensembl_ID) ->up2
merge(relationship,down,by.x="EnsemblId",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(EnsemblId %in% down$ensembl_ID) ->down2
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
dotplot(eReactome_up, showCategory=20)
enrichPathway(gene=down2$GeneID,organism="human",pvalueCutoff=0.05, readable=T)->eReactome_down
dotplot(eReactome_down, showCategory=20)->Reactome_down_plot

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



#heatmap for down-gene(p<0.01)----------------------------------------------------------
#dplyr::filter(down,p_value<=0.01) %>%
#  dplyr::filter(diff.avg<=-2)-> down_for_map
#rownames(down_for_map)=down_for_map[,3]
#heatmap(as.matrix(down_for_map[,4:(ncol(down_for_map)-4)]),Colv=NA,ColSideColors=c(rep("purple", 44), rep("orange", 117)),col=colorRampPalette(c("green", "black","red"))(256),cexRow = 0.3,cexCol = 0.2)

#Heatmap(ComplexHeatmap)
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/up.txt",sep="\t",header = T,as.is = TRUE) ->up
rownames(up)=up[,1]
up[order(up$p_value,decreasing=T),]->up2#diff.avg ordered from large to small
up2[,4:(ncol(up2)-4)] -> up2
up2[1:10,] %>% apply(2,function(x) scale(x)) ->a
rownames(a)=rownames(up2)[1:10]
Heatmap(as.matrix(a),cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2))


c(rep("1",44),rep("0",117)) %>% as.matrix() ->label
up3=rbind(label[,1],up2)
up3=as.data.frame(t(up3))
annot_up <- data.frame(response_VS_none = up3$`1`)
colors = list(response_VS_none = c("1" = "green", "0" = "gray"))
ha <- HeatmapAnnotation(annot_up, col = colors)
Heatmap(as.matrix(a),column_names_gp = gpar(fontsize = 2),top_annotation = ha,cluster_columns = FALSE)#Heatmap for ten top P

rownames(down2)=down2$EnsemblId
down2[,-c(1,2,3,(ncol(down2)-3),(ncol(down2)-2),(ncol(down2)-1),(ncol(down2)))]->a
annotation_col = data.frame(GeneClass = factor(rep(c("response", "non-response"), c(length(response_expression),length(non_response_expression)))))
rownames(annotation_col)=colnames(a)
pheatmap(a, annotation_col = annotation_col,cluster_cols = FALSE,scale="column")

