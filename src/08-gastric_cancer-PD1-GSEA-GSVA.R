library(dplyr)
library(magrittr)
library(fgsea)
library(reactome.db)


#fgsea in reactome------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/all_DEG.txt",header = T,as.is = TRUE) -> gastric_cancer_PD1_DEG
tibble::rownames_to_column(gastric_cancer_PD1_DEG) %>%
  dplyr::select(rowname,logFC)->EnsemblID_logFC
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,EnsemblID_logFC,by.x="Ensembl_ID",by.y="rowname")%>%
  dplyr::select(GeneID,logFC)->GeneID_logFC

gene_list=GeneID_logFC$logFC
names(gene_list)=GeneID_logFC$GeneID
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
write.table(sig_Reactome,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/GSEA/sig_Reactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#

# #plot the single significantly enriched pathway
# plotEnrichment(my_pathways[['Extracellular matrix organization']],
#                gene_list) + labs(title='Extracellular matrix organization')->ECM_organization
# ggsave(
#   filename = 'Down_ECM_organization.pdf',
#   plot = ECM_organization,
#   device = 'pdf',
#   path = '/data/liull/immune-checkpoint-blockade/GSEA/gastric_cancer_PD1/',
#   width = 12,
#   height = 8
# )
# 

#Plot the top 10 pathways enriched at the top and bottom of the ranked list, respectively
topPathwaysUp <- fgsea_reactome[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgsea_reactome[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(my_pathways[topPathways], ordered_gene_list, fgsea_reactome, 
              gseaParam = 0.5)->p

ggsave(
  filename = 'gastric_cancer_GSEA_top20.pdf',
  plot = p,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/GSEA/',
  width = 22,
  height = 8
)
rm(list = ls())


#GSVA-------------------------------------------------------------------
library(GSVA)
library(readxl)
library(magrittr)
library(GSEABase)
library(edgeR)
library(ComplexHeatmap)
library(circlize)

#load gene set,metadata and gene symbol,ensenml,entrezID relationship
file_path = "/data/liull/reference/GSEA-gmt"
genesets_c2<- getGmt(paste(file_path,"c2.all.v6.2.symbols.gmt",sep="/"))

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE,sep="\t") -> relationship

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="gastric cancer") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(Run,Response) ->gastric_cancer_PD1

dplyr::filter(gastric_cancer_PD1,Response %in% c("CR","PR"))$Run ->response_ids
dplyr::filter(gastric_cancer_PD1,Response %in% c("SD","PD"))$Run ->non_response_ids

#load gastric cancer expr and translate Ensembl id to symbol

read.table("/data/liull/immune-checkpoint-blockade/expression/all_count_expression_2.txt",header = T,as.is = TRUE) ->all_expression

all_expression %>%
  dplyr::select(gene_id,response_ids,non_response_ids)->gastric_cancer_expr
row.names(gastric_cancer_expr)=gastric_cancer_expr[,1]
gastric_cancer_expr=gastric_cancer_expr[,-1]

DGEList_expr <- DGEList(counts=gastric_cancer_expr)
normalized_expr <- calcNormFactors(DGEList_expr, method="upperquartile")
normalized_loggedCPM_expr = cpm(normalized_expr, log=TRUE, prior.count=2)

keep <- rowSums(normalized_loggedCPM_expr>0) >= 2
normalized_loggedCPM_expr <- normalized_loggedCPM_expr[keep,]#the same normalized(log2CPM) as melanoma

tibble::rownames_to_column(as.data.frame(normalized_loggedCPM_expr)) %>%
  dplyr::filter(rowname %in% relationship$Ensembl_ID) %>%
  merge(relationship,.,by.x="Ensembl_ID",by.y="rowname")%>%
  dplyr::select(Symbol,response_ids,non_response_ids)->expression

factors=factor(expression$Symbol)
merged_expression=tapply(expression[,2],factors,median)
for (i in 3:ncol(expression)) {
  temp=tapply(expression[,i],factors,median)
  merged_expression=cbind(merged_expression,temp)
}
colnames(merged_expression)=colnames(expression)[2:ncol(expression)]


#make GSVA score for melanoma_CTLA4
GSVA_score_c2 <- gsva(merged_expression, genesets_c2, min.sz=1, max.sz=999999, method="zscore",kcdf="Gaussian", abs.ranking=FALSE, verbose=TRUE)


group_list <- factor(c(rep("response",length(response_ids)), rep("non_response",length(non_response_ids))))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(GSVA_score_c2)

fit <- lmFit(GSVA_score_c2, design)
fit2 <- eBayes(fit)
output <- topTable(fit2, coef=2, n=Inf)
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC>1)->up#157
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC< -1)->down#32


#heatmap

rbind(up,down)->all_gene_sets
write.table(all_gene_sets,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/GSVA/c2_sig_sets.txt",sep="\t",quote=FALSE,col.names = TRUE,row.names = FALSE)

tibble::rownames_to_column(as.data.frame(GSVA_score_c2))%>%
  dplyr::filter(rowname %in% all_gene_sets$rowname)->GSVA_score_heatmap
rownames(GSVA_score_heatmap)=GSVA_score_heatmap$rowname
GSVA_score_heatmap=GSVA_score_heatmap[,-1]

apply(GSVA_score_heatmap, 1, scale) ->scaled_GSVA_score
rownames(scaled_GSVA_score)=colnames(GSVA_score_heatmap)
scaled_GSVA_score=t(scaled_GSVA_score)


df = data.frame(type = c(rep("response", length(response_ids)), rep("non_response", length(non_response_ids))))
ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))

pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/GSVA/heatmap.pdf")
Heatmap(scaled_GSVA_score,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 2),col=colorRamp2(c(-4, 0, 4), c("green", "black", "red")))
dev.off()
