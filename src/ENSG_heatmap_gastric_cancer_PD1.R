read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/all_DEG.txt",header = T,as.is = TRUE) ->all_DEG
tibble::rownames_to_column(all_DEG) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC>1)->up#668
tibble::rownames_to_column(all_DEG) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC< -1)->down#875

dplyr::filter(up,rowname %in% grep("ENSG",up$rowname,value=T))->up_ENSG
dplyr::filter(down,rowname %in% grep("ENSG",down$rowname,value = T)) ->down_ENSG


readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="gastric cancer") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

#select response and non_response's sample id and project id-------------------
dplyr::filter(metadata,Response %in% c("CR","PR")) ->response
dplyr::filter(metadata,Response %in% c("SD","PD")) ->non_response
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


#heatmap for ENSG
rbind(up_ENSG,down_ENSG)->all_genes
tibble::rownames_to_column(as.data.frame(normalized_loggedCPM_expr)) %>% 
  dplyr::filter(rowname %in% all_genes$rowname)->expr_heatmap
rownames(expr_heatmap)=expr_heatmap$rowname
expr_heatmap=expr_heatmap[,-1]

apply(expr_heatmap, 1, scale) ->scaled_expr
rownames(scaled_expr)=colnames(expr_heatmap)
scaled_expr=t(scaled_expr)


df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))

pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/heatmap_ENSG.pdf")
Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 0.5),col=colorRamp2(c(-4, 0, 4), c("green", "black", "red")))->origin_heatmap
dev.off()

#second heatmap
# > sum(rowSums(scaled_expr>3))
# [1] 89
# > sum(rowSums(scaled_expr< -3))
# [1] 242


new_scaled_expr <- scaled_expr[row_order(origin_heatmap)[[1]],]

for(i in 1:ncol(new_scaled_expr)) {
  m <- which(new_scaled_expr[,i]>3)
  new_scaled_expr[m,i] <- 3
  n <- which(new_scaled_expr[,i]<(-3))
  new_scaled_expr[n,i] <- (-3)
}
pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/heatmap_ENSG_2.pdf")
Heatmap(new_scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,cluster_rows = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 1),col=colorRamp2(c(-3, 0, 3), c("green", "black", "red")))
dev.off()

#heatmap for NONHSAG
dplyr::filter(up,rowname %in% grep("NONHSAG",up$rowname,value=T))->up_NONHSAG
dplyr::filter(down,rowname %in% grep("NONHSAG",down$rowname,value = T)) ->down_NONHSAG
rbind(up_NONHSAG,down_NONHSAG)->NONHSAG_genes
tibble::rownames_to_column(as.data.frame(normalized_loggedCPM_expr)) %>% 
  dplyr::filter(rowname %in% NONHSAG_genes$rowname)->NONHSAG_expr_heatmap
rownames(NONHSAG_expr_heatmap)=NONHSAG_expr_heatmap$rowname
NONHSAG_expr_heatmap=NONHSAG_expr_heatmap[,-1]

apply(NONHSAG_expr_heatmap, 1, scale) ->NONHSAG_scaled_expr
rownames(NONHSAG_scaled_expr)=colnames(NONHSAG_expr_heatmap)
NONHSAG_scaled_expr=t(NONHSAG_scaled_expr)


df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))

pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/heatmap_NONHSAG.pdf")
Heatmap(NONHSAG_scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 0.5),col=colorRamp2(c(-4, 0, 4), c("green", "black", "red")))
dev.off()





# read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
# 
# merge(relationship,up,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
#   dplyr::filter(Ensembl_ID %in% up$rowname) ->up2
# up2<- up2[order(up2$logFC,decreasing = TRUE),]
# merge(relationship,down,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
#   dplyr::filter(Ensembl_ID %in% down$rowname) ->down2
# down2<- down2[order(down2$logFC),]
# write.table(up2,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_symbol_ordered.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)
# write.table(down2,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_symbol_ordered.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)
