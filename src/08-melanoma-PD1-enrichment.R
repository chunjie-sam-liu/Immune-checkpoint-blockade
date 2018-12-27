library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)


read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/all_values.txt",sep="\t",header = T,as.is = TRUE) ->all_exp_data
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/up_gene_list.txt",sep="\t",header = T,as.is = TRUE) ->up_list
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/down_gene_list.txt",sep="\t",header = T,as.is = TRUE) ->down_list


trans_up_list=bitr(up_list$Gene.name,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
up_go_ALL <- enrichGO(gene = trans_up_list$ENTREZID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE)
#write.csv(summary(up_go_ALL),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/up_enrich.csv",row.names =FALSE)
up=as.data.frame(up_go_ALL)
GO_up_gene_list=unique(sapply(up$geneID, function(x) strsplit(x,"[/]")))

trans_down_list=bitr(down_list$Gene.name,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
down_go_ALL <- enrichGO(gene = trans_down_list$ENTREZID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE)
#write.csv(summary(down_go_ALL),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/down_enrich.csv",row.names =FALSE)
down=as.data.frame(down_go_ALL)
GO_down_gene_list=unique(sapply(down$geneID, function(x) strsplit(x,"[/]")))


expression=all_exp_data[,-((ncol(all_exp_data)-3):ncol(all_exp_data))]
dplyr::mutate(expression,symbol=bitr(expression$ensembl_ID,fromType="ENSEMBL",toType="SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL) -> expression_2

rownames(expression)=expression[,1]
expression=expression[,-1]
es.max <- gsva(expression, up_gene, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
