library(magrittr)
library(readr)
library(readxl)
library(dplyr)
library(sva)

#筛选出gastric cancer的RNA-seq的anti-PD1的数据的应答信息-------------
readxl::read_excel("/data/liull/immune-checkpoint-blockade/04-all-metadata.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="gastric cancer") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

#select response and non_response's sample id and project id-------------------
rbind(dplyr::filter(metadata,Response=="CR"),dplyr::filter(metadata,Response=="PR")) ->response
rbind(dplyr::filter(metadata,Response=="SD"),dplyr::filter(metadata,Response=="PD")) ->non_response

#make difference----------------------------
read.table("/data/liull/immune-checkpoint-blockade/expression/all_expression.txt",sep="\t",header = T,as.is = TRUE) ->data
all_expression=dplyr::select(data,gene_id,response$Run,non_response$Run)
#order the expression profile by response and nonresponse to t.test


######
Sum_NA=apply(all_expression,1,function(x) sum(is.na(x)))
NA_IDs=which(Sum_NA>=(length(all_expression)/4))
all_expression2=all_expression[-NA_IDs,]#delete the gene has more than 1/4 samples' NA

row.names(all_expression2)=all_expression2[,1]
all_expression2=all_expression2[,-1]
dim(all_expression2)#make rownames to avoid of sum wrong

Sum_zero=apply(all_expression2,1,function(x) sum(as.numeric(x),na.rm = TRUE))
IDs_zero=which(Sum_zero==0)
all_expression3=all_expression2[-IDs_zero,]#delete the gene has all 0.000 depression

for(i in 1:length(all_expression3)) {
  all_expression3[is.na(all_expression3[, i]), i] <- mean(all_expression3[, i], na.rm = T)
}#replace NA to mean of its sample expression
######

#t.test
avg.R=apply(all_expression3,1,function(x) mean(x[1:nrow(response)]))
avg.NR=apply(all_expression3,1,function(x) mean(x[(nrow(response)+1):length(all_expression3)]))
FC=apply(all_expression3,1,function(x) (mean(x[1:nrow(response)])+0.01)/(mean(x[(nrow(response)+1):length(all_expression3)])+0.01))
p_value=apply(all_expression3,1,function(x) t.test(x[1:nrow(response)],x[(nrow(response)+1):length(all_expression3)])$p.value)
result=as.data.frame(cbind(all_expression3,avg.R,avg.NR,FC,p_value))

result=cbind(rownames(result),result)
rownames(result)=NULL
colnames(result)[1]="ensembl_ID"
result$ensembl_ID <- as.character(result$ensembl_ID)

dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
  dplyr::filter(FC>=2) -> up
dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
  dplyr::filter(FC<=0.5) -> down

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,up,by.x="EnsemblId",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(EnsemblId %in% up$ensembl_ID) ->up2
merge(relationship,down,by.x="EnsemblId",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(EnsemblId %in% down$ensembl_ID) ->down2


write.table(up2,"/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/up.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(down2,"/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/down.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)


#GO enrichment
enrichGO(gene = up2$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,readable = TRUE) %>% as.data.frame() ->up_enrichGO#dim(up_enrichGO) 3
enrichGO(gene = down2$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,readable = TRUE) %>% as.data.frame() ->down_enrichGO#55
write.table(up_enrichGO,"/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/up_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(down_enrichGO,"/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/down_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#draw a map for most 10 significant enriched
down_enrichGO=down_enrichGO[order(down_enrichGO$p.adjust),]
go_enrich_df<- data.frame(GO=down_enrichGO$Description,p.adjust=down_enrichGO$p.adjust)[1:10,]

ggplot(data=go_enrich_df, aes(x=GO, y=-log10(p.adjust))) +
  geom_bar(stat="identity", width=0.6) + coord_flip()+theme_bw()+
  labs(title = "The Most Enriched GO Terms") -> GO_enrich_pdf
ggsave(
  filename = 'gastric_cancer_PD1_GO_enrich.pdf',
  plot = GO_enrich_pdf,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer',
  width = 8,
  height = 6.8
)

#KEGG enrichment
enrichKEGG(gene=up2$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH") %>%
  as.data.frame()->up_enrichKEGG#0
enrichKEGG(gene=down2$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH") %>%
  as.data.frame()->down_enrichKEGG
write.table(down_enrichKEGG,"/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/down_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
