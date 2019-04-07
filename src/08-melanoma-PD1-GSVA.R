library(GSVA)
library(readxl)
library(magrittr)
library(GSEABase)

#load gene set,metadata and gene symbol,ensenml,entrezID relationship---------------------------------------
file_path = "/data/liull/reference/GSEA-gmt"
genesets_c2<- getGmt(paste(file_path,"c2.all.v6.2.symbols.gmt",sep="/"))

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE,sep="\t") -> relationship

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(Run,Response) ->melanoma_PD1
dplyr::filter(melanoma_PD1,Response %in% c("CR","PR","R","PRCR"))$Run ->response_ids
dplyr::filter(melanoma_PD1,Response %in% c("SD","PD","NR"))$Run ->non_response_ids

#load melanoma_PD1_removed_batch_expression and translate Ensembl id to symbol--------------------------------------------

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/PD1_removed_batch_expression.txt",header = T,as.is = TRUE) ->all_expression

tibble::rownames_to_column(all_expression) %>%
  dplyr::filter(rowname %in% relationship$Ensembl_ID) %>%
  merge(relationship[,2:3],.,by.x="Ensembl_ID",by.y="rowname")%>%
  dplyr::select(Symbol,response_ids,non_response_ids)->expression
rownames(expression) <- expression[,1]
expression <- expression[,-1] 



#make GSVA score for melanoma_PD1
GSVA_score_c2 <- gsva(data.matrix(expression), genesets_c2, min.sz=1, max.sz=999999, method="zscore",kcdf="Gaussian", abs.ranking=FALSE, verbose=TRUE)


group_list <- factor(c(rep("response",length(response_ids)), rep("non_response",length(non_response_ids))))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(GSVA_score_c2)

fit <- lmFit(GSVA_score_c2, design)
fit2 <- eBayes(fit)
output <- topTable(fit2, coef=2, n=Inf)
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC>1)->up
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC< -1)->down


#heatmap

rbind(up,down)->all_gene_sets
write.table(all_gene_sets,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/GSVA/c2_sig_sets.txt",sep="\t",quote=FALSE,col.names = TRUE,row.names = FALSE)

tibble::rownames_to_column(as.data.frame(GSVA_score_c2))%>%
  dplyr::filter(rowname %in% all_gene_sets$rowname)->GSVA_score_heatmap
rownames(GSVA_score_heatmap)=GSVA_score_heatmap$rowname
GSVA_score_heatmap=GSVA_score_heatmap[,-1]

apply(GSVA_score_heatmap, 1, scale) ->scaled_GSVA_score
rownames(scaled_GSVA_score)=colnames(GSVA_score_heatmap)
scaled_GSVA_score=t(scaled_GSVA_score)


df = data.frame(type = c(rep("response", length(response_ids)), rep("non_response", length(non_response_ids))))
ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))

pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/GSVA/heatmap.pdf")
Heatmap(scaled_GSVA_score,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 1),col=colorRamp2(c(-4, 0, 4), c("green", "black", "red")))
dev.off()