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

GO_up_gene_list=unique(sapply(up_go_ALL$geneID, function(x) strsplit(x,"[/]")))
GO_up_gene_list_2=GO_up_gene_list
for (i in 1:nrow(trans_up_list)) {
  for (j in 1:length(GO_up_gene_list)) {
    gsub(trans_up_list$Gene.name[i],trans_up_list$Gene.stable.ID[i],GO_up_gene_list_2[j]) ->GO_up_gene_list_2[j]
  }
}

GO_up_gene_list_2=as.character(GO_up_gene_list_2)
GO_up_gene_list_3=list()
length(GO_up_gene_list_3)=length(GO_up_gene_list_2)
for (i in 1:length(GO_up_gene_list_2)) {
  number=length(gregexpr("ENSG",GO_up_gene_list_2[i])[[1]])
  #print(number)
  for(j in 1:number){
    n=gregexpr("ENSG",GO_up_gene_list_2[i])[[1]][j]
    #print(n)
    GO_up_gene_list_3[[i]][j]=substr(GO_up_gene_list_2[i], n, n+14)
  }
}#get all of the gene set as list form(no enrichment cutoff)


#down regulated gene GO enrichment
bitr(down_list$Gene.name,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db") %>%
  merge(down_list,.,by.x="Gene.name",by.y="SYMBOL",all=T) ->trans_down_list #Genename  GeneEnsemblID  ENTREZID

enrichGO(gene = trans_down_list$ENTREZID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE) %>%
  as.data.frame() -> down_go_ALL


GO_down_gene_list=unique(sapply(down_go_ALL$geneID, function(x) strsplit(x,"[/]")))
GO_down_gene_list_2=GO_down_gene_list
for (i in 1:nrow(trans_down_list)) {
  for (j in 1:length(GO_down_gene_list)) {
    gsub(trans_down_list$Gene.name[i],trans_down_list$Gene.stable.ID[i],GO_down_gene_list_2[j]) ->GO_down_gene_list_2[j]
  }
}
GO_down_gene_list_2=as.character(GO_down_gene_list_2)
GO_down_gene_list_3=list()
length(GO_down_gene_list_3)=length(GO_down_gene_list_2)
for (i in 1:length(GO_down_gene_list_2)) {
  number=length(gregexpr("ENSG",GO_down_gene_list_2[i])[[1]])
  #print(number)
  for(j in 1:number){
    n=gregexpr("ENSG",GO_down_gene_list_2[i])[[1]][j]
    #print(n)
    GO_down_gene_list_3[[i]][j]=substr(GO_down_gene_list_2[i], n, n+14)
  }
}#get all of the gene set(no enrichment cutoff)


expression=all_exp_data[,-((ncol(all_exp_data)-3):ncol(all_exp_data))]
rownames(expression)=expression[,1]
expression=expression[,-1]

up_gsva <- gsva(as.matrix(expression), GO_up_gene_list_3, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
rownames(up_gsva)=GO_up_gene_list_3
down_gsva <- gsva(as.matrix(expression), GO_down_gene_list_3, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
rownames(down_gsva)=GO_down_gene_list_3

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

dplyr::select(as.data.frame(down_gsva),response_ids) ->down_response_gsva
dplyr::select(as.data.frame(down_gsva),non_response_ids) ->down_non_response_gsva
ordered_down_gsva=cbind(down_response_gsva,down_non_response_gsva)


