readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(SRA_Study,Run,Response,Biopsy_Time) ->metadata

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/PD1_removed_batch_expression.txt",header = T,as.is = TRUE) ->combat_edata

#DEG by limma
dplyr::filter(metadata,Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(Response %in% c("CR","PR","PRCR","R"))-> response
dplyr::filter(metadata,Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(Response %in% c("SD","PD","NR")) -> non_response

dplyr::select(as.data.frame(combat_edata),response$Run,non_response$Run)->ordered_combat_edata

keep <- rowSums(ordered_combat_edata>0) >= 2
ordered_combat_edata <- ordered_combat_edata[keep,]
#delete the gene has less than 2 sample exression CPM<1(log2CPM<0)

group_list <- factor(c(rep("response",nrow(response)), rep("non_response",nrow(non_response))))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(ordered_combat_edata)

fit <- lmFit(ordered_combat_edata, design)
fit2 <- eBayes(fit)
output <- topTable(fit2, coef=2, n=Inf)
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC>1)->up
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC< -1)->down

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship

merge(relationship,up,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% up$rowname) ->up2
up2<- up2[order(up2$logFC,decreasing = TRUE),]
merge(relationship,down,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% down$rowname) ->down2
down2<- down2[order(down2$logFC),]
write.table(up2,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_up_symbol.txt",quote = FALSE,row.names = TRUE,col.names = TRUE)
write.table(down2,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_down_symbol.txt",quote = FALSE,row.names = TRUE,col.names = TRUE)


dplyr::filter(up,rowname %in% grep("ENSG",up$rowname,value=T))->up_ENSG
dplyr::filter(down,rowname %in% grep("ENSG",down$rowname,value = T)) ->down_ENSG

rbind(up_ENSG,down_ENSG)->all_genes
tibble::rownames_to_column(ordered_combat_edata) %>% 
  dplyr::filter(rowname %in% all_genes$rowname)->expr_heatmap
rownames(expr_heatmap)=expr_heatmap$rowname
expr_heatmap=expr_heatmap[,-1]

apply(expr_heatmap, 1, scale) ->scaled_expr
rownames(scaled_expr)=colnames(expr_heatmap)
scaled_expr=t(scaled_expr)


df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))

pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/heatmap_ENSG.pdf")
Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 1),col=colorRamp2(c(-4, 0, 4), c("green", "black", "red")))
dev.off()

# second heatmap

new_scaled_expr <- scaled_expr[row_order(origin_heatmap)[[1]],]

for(i in 1:ncol(new_scaled_expr)) {
  m <- which(new_scaled_expr[,i]>2)
  new_scaled_expr[m,i] <- 2
  n <- which(new_scaled_expr[,i]<(-2))
  new_scaled_expr[n,i] <- (-2)
}
pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/heatmap_ENSG_2.pdf")
Heatmap(new_scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,cluster_rows = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 1),col=colorRamp2(c(-2, 0, 2), c("green", "black", "red")))
dev.off()



#filtered at first
tibble::rownames_to_column(combat_edata) %>%
  dplyr::filter(rowname %in% grep("ENSG",rownames(combat_edata),value=T))->ENSG_combat_edata
dplyr::select(ENSG_combat_edata,rowname,response$Run,non_response$Run)->ordered_ENSG_combat_edata
rownames(ordered_ENSG_combat_edata)=ordered_ENSG_combat_edata$rowname
ordered_ENSG_combat_edata=ordered_ENSG_combat_edata[,-1]

keep <- rowSums(ordered_ENSG_combat_edata>0) >= 2
ordered_ENSG_combat_edata <- ordered_ENSG_combat_edata[keep,]
#delete the gene has less than 2 sample exression CPM<1(log2CPM<0)

group_list <- factor(c(rep("response",nrow(response)), rep("non_response",nrow(non_response))))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(ordered_ENSG_combat_edata)

fit <- lmFit(ordered_ENSG_combat_edata, design)
fit2 <- eBayes(fit)
output <- topTable(fit2, coef=2, n=Inf)
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC>1)->up
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC< -1)->down

              