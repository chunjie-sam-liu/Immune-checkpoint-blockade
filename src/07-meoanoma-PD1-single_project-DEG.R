library(ComplexHeatmap)
library(circlize)

#metadata-----------------------------------------------------------------------------------------------------------------
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response)%>%
  dplyr::filter(Response != "NE")->metadata


#expr--------------------------------------------------------------------------------------------------------------

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_count_expr.txt",
           header = T,as.is = TRUE) %>%
  dplyr::select(metadata$Run)->PD1_count_expr

#DEG for each single project-----------------------------------------------------------------------------------
c("SRP070710","SRP150548","SRP094781")->project_id
for (i in 1:length(project_id)) {
  
  dplyr::filter(metadata,SRA_Study == project_id[i])->SRP_metadata
  
  SRP_metadata%>%
    dplyr::filter(Response %in% c("CR","PR","R"))-> response
  SRP_metadata%>%
    dplyr::filter(Response %in% c("SD","PD","NR")) -> non_response
  
  dplyr::select(PD1_count_expr,response$Run,non_response$Run)->ordered_expr
  
  keep <- rowSums(ordered_expr>0) >= 2
  ordered_expr <- ordered_expr[keep,]
  
  group_list <- factor(c(rep("response",nrow(response)), rep("non_response",nrow(non_response))))
  design <- model.matrix(~group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(ordered_expr)
  
  fit <- lmFit(ordered_expr, design)
  fit2 <- eBayes(fit)
  output <- topTable(fit2, coef=2, n=Inf)
  tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC>1) %>% tibble::column_to_rownames()->up
  tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC< -1)%>% tibble::column_to_rownames()->down
  
  write.table(output,paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],"all_DEG.txt",sep = "/"),
              quote = FALSE,sep="\t",row.names = TRUE,col.names = TRUE)
  write.table(up,paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],"up.txt",sep = "/"),
              quote = FALSE,row.names = TRUE,col.names = TRUE,sep = "\t")
  write.table(down,paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],"down.txt",sep = "/"),
              quote = FALSE,row.names = TRUE,col.names = TRUE,sep = "\t")

  
  
  c(rownames(up),rownames(down))->all_genes
  tibble::rownames_to_column(ordered_expr) %>% 
    dplyr::filter(rowname %in% all_genes)%>%
    tibble::column_to_rownames()->expr_heatmap
  
  apply(expr_heatmap, 1, scale) ->scaled_expr
  rownames(scaled_expr)=colnames(expr_heatmap)
  scaled_expr=t(scaled_expr)
  
  
  df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
  ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))
  
  pdf(file=paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],"heatmap.pdf",sep = "/"))
  Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 8),row_names_gp = gpar(fontsize = 0.5),col=colorRamp2(c(-4, 0, 4), c("green", "black", "red")))
  dev.off()
  
  
}


#DEG interaction------------------------------------------------------------------------------------------------
#up
read.table("/data/liull/immune-checkpoint-blockade/single_project_result/SRP070710/up.txt",header = T,as.is = TRUE)->SRP070710_up
read.table("/data/liull/immune-checkpoint-blockade/single_project_result/SRP094781/up.txt",header = T,as.is = TRUE)->SRP094781_up
read.table("/data/liull/immune-checkpoint-blockade/single_project_result/SRP150548/up.txt",header = T,as.is = TRUE)->SRP150548_up

intersect(rownames(SRP070710_up),rownames(SRP094781_up))
intersect(rownames(SRP070710_up),rownames(SRP150548_up))
intersect(rownames(SRP094781_up),rownames(SRP150548_up))

venn.diagram(list(SRP070710=rownames(SRP070710_up), SRP094781=rownames(SRP094781_up),SRP150548=rownames(SRP150548_up)),
                  filename=NULL,fill = c("cornflowerblue", "green", "yellow"))->venn_up
grid.draw(venn_up)

#down
read.table("/data/liull/immune-checkpoint-blockade/single_project_result/SRP070710/down.txt",header = T,as.is = TRUE)->SRP070710_down
read.table("/data/liull/immune-checkpoint-blockade/single_project_result/SRP094781/down.txt",header = T,as.is = TRUE)->SRP094781_down
read.table("/data/liull/immune-checkpoint-blockade/single_project_result/SRP150548/down.txt",header = T,as.is = TRUE)->SRP150548_down

intersect(rownames(SRP070710_down),rownames(SRP094781_down))
intersect(rownames(SRP070710_down),rownames(SRP150548_down))
intersect(rownames(SRP094781_down),rownames(SRP150548_down))# 0


venn.diagram(list(SRP070710=rownames(SRP070710_down), SRP094781=rownames(SRP094781_down),SRP150548=rownames(SRP150548_down)),
             filename=NULL,fill = c("cornflowerblue", "green", "yellow"))->venn_down
grid.draw(venn_down)


#GO,KEGG enrichment-------------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)


c("SRP070710","SRP150548","SRP094781")->project_id
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,sep="\t",as.is = TRUE) ->relationship

for (i in 1:length(project_id)) {
  
  #up
  read.table(paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],"up.txt",sep = "/"),header = T,as.is = TRUE)->up
  tibble::rownames_to_column(up)%>%
    merge(relationship,by.x="rowname",by.y="Symbol")->up
  
  ##GO
  enrichGO(gene = up$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "fdr",pvalueCutoff = 0.05,readable = TRUE)->ego_up
  write.table(as.data.frame(ego_up),paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],"GOenrich_up.txt",sep = "/"),
              quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

  DOSE::dotplot(ego_up, split="ONTOLOGY",showCategory=10) + facet_grid(ONTOLOGY~., scale="free")->ego_up_plot
  ggsave(
    filename = 'GOenrich_up.pdf',
    plot = ego_up_plot,
    device = 'pdf',
    path = paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],sep = "/"),
    width = 12,
    height = 8
  )
  
  ##KEGG
  enrichKEGG(gene=up$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "fdr") ->ekegg_up
  write.table(as.data.frame(ekegg_up),paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],"KEGGenrich_up.txt",sep = "/"),
              quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
  
  DOSE::dotplot(ekegg_up, showCategory=10)->KEGG_up_plot
  ggsave(
    filename = 'KEGGenrich_up.pdf',
    plot = KEGG_up_plot,
    device = 'pdf',
    path = paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],sep = "/"),
    width = 12,
    height = 8
  )
  
  
  #down
  read.table(paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],"down.txt",sep = "/"),header = T,as.is = TRUE)->down
  tibble::rownames_to_column(down)%>%
    merge(relationship,by.x="rowname",by.y="Symbol")->down
  
  ##GO
  enrichGO(gene = down$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "fdr",pvalueCutoff = 0.05,readable = TRUE)->ego_down
  write.table(as.data.frame(ego_down),paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],"GOenrich_down.txt",sep = "/"),
              quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
  
  DOSE::dotplot(ego_down, split="ONTOLOGY",showCategory=10) + facet_grid(ONTOLOGY~., scale="free")->ego_down_plot
  ggsave(
    filename = 'GOenrich_down.pdf',
    plot = ego_down_plot,
    device = 'pdf',
    path = paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],sep = "/"),
    width = 12,
    height = 8
  )
  
  #KEGG
  enrichKEGG(gene=down$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "fdr") ->ekegg_down
  write.table(as.data.frame(ekegg_down),paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],"KEGGenrich_down.txt",sep = "/"),
              quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
  
  DOSE::dotplot(ekegg_down, showCategory=10)->KEGG_down_plot
  ggsave(
    filename = 'KEGGenrich_down.pdf',
    plot = KEGG_down_plot,
    device = 'pdf',
    path = paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],sep = "/"),
    width = 12,
    height = 8
  )
  
}


#GSEA----------------------------------------------------------------------------
library(fgsea)
library(reactome.db)

c("SRP070710","SRP150548","SRP094781")->project_id
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,sep="\t",as.is = TRUE) ->relationship

for (i in 1:length(project_id)) {
  
  read.table(paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],"all_DEG.txt",sep = "/"),header = T,as.is = TRUE)->all_genes
  tibble::rownames_to_column(all_genes)%>%
    merge(relationship,by.x="rowname",by.y="Symbol")->all_genes
  
  gene_list=all_genes$logFC
  names(gene_list)=all_genes$GeneID
  ordered_gene_list <- gene_list[order(gene_list)]
  
  
  my_pathways <- reactomePathways(names(ordered_gene_list))
  fgsea_reactome <- fgsea(pathways = my_pathways, 
                          stats = ordered_gene_list,
                          minSize=15,
                          maxSize=500,
                          nperm=100000)
  head(fgsea_reactome[order(pval), ])
  sum(fgsea_reactome[, padj < 0.05])
  
  dplyr::filter(fgsea_reactome,padj< 0.05)->sig_Reactome
  data.frame(lapply(sig_Reactome,as.character), stringsAsFactors=FALSE)->sig_Reactome
  write.table(sig_Reactome,paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],"GSEA.txt",sep = "/"),
              quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
  
  
  topPathwaysUp <- fgsea_reactome[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgsea_reactome[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(my_pathways[topPathways], ordered_gene_list, fgsea_reactome, 
                gseaParam = 0.5)->p
  
  ggsave(
    filename = 'GSEA_top20.pdf',
    plot = p,
    device = 'pdf',
    path = paste("/data/liull/immune-checkpoint-blockade/single_project_result",project_id[i],sep = "/"),
    width = 22,
    height = 8
  )
  
}




