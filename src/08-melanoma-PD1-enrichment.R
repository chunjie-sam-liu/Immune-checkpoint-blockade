library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(readxl)

read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/all_values.txt",sep="\t",header = T,as.is = TRUE) ->all_exp_data
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/up_gene_list.txt",sep="\t",header = T,as.is = TRUE) ->up_list
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/down_gene_list.txt",sep="\t",header = T,as.is = TRUE) ->down_list


#up regulated gene GO enrichment
bitr(up_list$Gene.name,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db") %>%
  merge(up_list,.,by.x="Gene.name",by.y="SYMBOL",all=T) ->trans_up_list #Genename  GeneEnsemblID  ENTREZID


enrichGO(gene = trans_up_list$ENTREZID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE) %>%
  as.data.frame() -> up_go_ALL

GO_up_gene_list=sapply(up_go_ALL$geneID, function(x) strsplit(x,"[/]"))
#get all of the gene set as list form(no enrichment cutoff)


#down regulated gene GO enrichment
bitr(down_list$Gene.name,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db") %>%
  merge(down_list,.,by.x="Gene.name",by.y="SYMBOL",all=T) ->trans_down_list #Genename  GeneEnsemblID  ENTREZID

enrichGO(gene = trans_down_list$ENTREZID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE) %>%
  as.data.frame() -> down_go_ALL


GO_down_gene_list=sapply(down_go_ALL$geneID, function(x) strsplit(x,"[/]"))
#get all of the gene set(no enrichment cutoff)


expression=all_exp_data[,-((ncol(all_exp_data)-3):ncol(all_exp_data))]

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE,sep="\t") -> relationship
dplyr::filter(expression,ensembl_ID %in% relationship$EnsemblId) ->expression_1
merge(relationship[,2:3],expression_1,by.x="EnsemblId",by.y="ensembl_ID") -> expression_2

rownames(expression_2)=expression_2[,2]
expression_2=expression_2[,-c(1,2)]


gsva(as.matrix(expression_2), GO_up_gene_list, mx.diff=FALSE, verbose=FALSE, parallel.sz=1) %>% 
  unique() -> up_gsva
gsva(as.matrix(expression_2), GO_down_gene_list, mx.diff=FALSE, verbose=FALSE, parallel.sz=1) %>%
  unique() -> down_gsva

readxl::read_excel("/data/liull/immune-checkpoint-blockade/04-all-metadata.xlsx",col_names = TRUE,sheet="SRA") ->metadata 

metadata %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(Run,Response) ->melanoma_PD1

#select response and non_response's sample id and project id
dplyr::filter(melanoma_PD1,Response %in% c("CR","PR","R"))$Run ->response_ids
dplyr::filter(melanoma_PD1,Response %in% c("SD","PD","NR"))$Run ->non_response_ids
#order gsva result by response and non-response
dplyr::select(as.data.frame(up_gsva),response_ids) ->up_response_gsva
dplyr::select(as.data.frame(up_gsva),non_response_ids) ->up_non_response_gsva
ordered_up_gsva=cbind(up_response_gsva,up_non_response_gsva)

response_Mean_up=apply(ordered_up_gsva,1,function(x) mean(x[1:ncol(up_response_gsva)]))
non_response_Mean_up=apply(ordered_up_gsva,1,function(x) mean(x[(ncol(up_response_gsva)+1):length(ordered_up_gsva)]))
FC_up=apply(ordered_up_gsva,1,function(x) (mean(x[1:ncol(up_response_gsva)])+0.01)/(mean(x[(ncol(up_response_gsva)+1):length(ordered_up_gsva)])+0.01))
p_value_up=apply(ordered_up_gsva,1,function(x) t.test(x[1:ncol(up_response_gsva)],x[(ncol(up_response_gsva)+1):length(ordered_up_gsva)])$p.value)
result_up=as.data.frame(cbind(ordered_up_gsva,response_Mean_up,non_response_Mean_up,FC_up,p_value_up))


dplyr::select(as.data.frame(down_gsva),response_ids) ->down_response_gsva
dplyr::select(as.data.frame(down_gsva),non_response_ids) ->down_non_response_gsva
ordered_down_gsva=cbind(down_response_gsva,down_non_response_gsva)

response_Mean_down=apply(ordered_down_gsva,1,function(x) mean(x[1:ncol(down_response_gsva)]))
non_response_Mean_down=apply(ordered_down_gsva,1,function(x) mean(x[(ncol(down_response_gsva)+1):length(ordered_down_gsva)]))
FC_down=apply(ordered_down_gsva,1,function(x) (mean(x[1:ncol(down_response_gsva)])+0.01)/(mean(x[(ncol(down_response_gsva)+1):length(ordered_down_gsva)])+0.01))
p_value_down=apply(ordered_down_gsva,1,function(x) t.test(x[1:ncol(down_response_gsva)],x[(ncol(down_response_gsva)+1):length(ordered_down_gsva)])$p.value)
result_down=as.data.frame(cbind(ordered_down_gsva,response_Mean_down,non_response_Mean_down,FC_down,p_value_down))
