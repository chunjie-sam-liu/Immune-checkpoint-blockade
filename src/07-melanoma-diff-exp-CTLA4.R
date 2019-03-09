library(magrittr)
library(readr)
library(readxl)
library(dplyr)

#load melanoma_CTLA4_removed_batch profile
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_CTLA4_removed_batch_expression.txt",header = T,as.is = TRUE)->combat_edata

#filter melanoma RNA-seq anti-CTLA4 metadata
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP")->dbGAP
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA")->SRA
rbind(dbGAP,SRA)%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

#select response and non_response's sample id and project id
dplyr::filter(metadata,Response %in% c("CR","PR","PRCR","R")) -> response
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response
dplyr::select(combat_edata,response$Run,non_response$Run)->ordered_melanoma_CTLA4

#orderthe expression profile by response and nonresponse to t.test

avg.R=apply(ordered_melanoma_CTLA4,1,function(x) median(x[1:nrow(response)]))
avg.NR=apply(ordered_melanoma_CTLA4,1,function(x) median(x[(nrow(response)+1):length(ordered_melanoma_CTLA4)]))
diff.avg=apply(ordered_melanoma_CTLA4,1,function(x) (median(x[1:nrow(response)])-median(x[(nrow(response)+1):length(ordered_melanoma_CTLA4)])))
p_value=apply(ordered_melanoma_CTLA4,1,function(x) t.test(x[1:nrow(response)],x[(nrow(response)+1):length(ordered_melanoma_CTLA4)])$p.value)
t_statistic=apply(ordered_melanoma_CTLA4,1,function(x) t.test(x[1:nrow(response)],x[(nrow(response)+1):length(ordered_melanoma_CTLA4)])$statistic)
result=as.data.frame(cbind(ordered_melanoma_CTLA4,avg.R,avg.NR,diff.avg,p_value,t_statistic))

result=cbind(rownames(result),result)
rownames(result)=NULL
colnames(result)[1]="ensembl_ID"
write.table(result,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/melanoma_CTLA4_DEG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#write all genes' difference

#filter significant different gene-------------------------------------
dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
  dplyr::filter(diff.avg>=2) -> up#53

dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
  dplyr::filter(diff.avg<=-2) -> down#98

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,up,by.x="EnsemblId",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(EnsemblId %in% up$ensembl_ID) ->up2
merge(relationship,down,by.x="EnsemblId",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(EnsemblId %in% down$ensembl_ID) ->down2
write.table(up2,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/up.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(down2,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/down.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#GO enrich for down-gene---------------------------------------------------------------------------------------
enrichGO(gene = up2$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.1,readable = TRUE)  ->ego_up
dotplot(ego_up)#0
enrichGO(gene = down2$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.1,readable = TRUE)->ego_down#
dotplot(ego_down,showCategory=20)->ego_down_plot
ggsave(
  filename = 'melanoma_CTLA4_down_GOenrich.pdf',
  plot = ego_down_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4',
  width = 12,
  height = 8
)

#write.table(as.data.frame(ego_up),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/up_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#0
write.table(as.data.frame(ego_down),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/down_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#1

#KEGG enrichment-----------------------------------------------------------------------
enrichKEGG(gene=up2$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH") ->ekegg_up#0
enrichKEGG(gene=down2$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH") ->ekegg_down#2
browseKEGG(ekegg_down, 'hsa04110')

write.table(down_enrichKEGG,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/down_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#Reactome enrichment------------------------------------------------------------
enrichPathway(gene=up2$GeneID,pvalueCutoff=0.1, readable=T)->eReactome_up#0
dotplot(eReactome_up)#0
enrichPathway(gene=down2$GeneID,pvalueCutoff=0.1, readable=T)->eReactome_down#5
dotplot(eReactome_down, showCategory=20)->Reactome_down_plot#3
ggsave(
  filename = 'melanoma_CTLA4_down_Reactome.pdf',
  plot = Reactome_down_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4',
  width = 12,
  height = 8
)

write.table(as.data.frame(eReactome_down),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/down_enrichReactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

