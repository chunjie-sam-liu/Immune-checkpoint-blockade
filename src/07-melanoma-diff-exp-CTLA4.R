library(magrittr)
library(readr)
library(readxl)
library(dplyr)

#筛选出黑色素瘤的RNA-seq的anti-CTLA4的数据的应答信息
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP")->dbGAP
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA")->SRA
rbind(dbGAP,SRA)%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

#select response and non_response's sample id and project id
dplyr::filter(metadata,Response %in% c("CR","PR","PRCR","R")) -> response
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response

#expression prepare for batch effect----------------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/expression/all_FPKM_expression_2.txt",sep="\t",header = T,as.is = TRUE) ->data1
Project=unique(metadata$SRA_Study)
dplyr::filter(metadata,SRA_Study==Project[1]) %>%
  dplyr::select(Run)  %>%
  as.matrix() %>%
  as.character()->Project1_id

dplyr::filter(metadata,SRA_Study==Project[2]) %>%
  dplyr::select(Run)  %>%
  as.matrix() %>%
  as.character()->Project2_id
expression1=dplyr::select(data1,gene_id,Project1_id,Project2_id)
#order the expression profile by project to remove batch effect

Sum_NA=apply(expression1,1,function(x) sum(is.na(x)))
NA_IDs=which(Sum_NA>=(length(expression1)/4))
expression2=expression1[-NA_IDs,]#delete the gene has more than 1/4 samples' NA

row.names(expression2)=expression2[,1]
expression2=expression2[,-1]
dim(expression2)#make rownames to avoid of sum wrong

Sum_zero=apply(expression2,1,function(x) sum(as.numeric(x),na.rm = TRUE))
IDs_zero=which(Sum_zero==0)
expression3=expression2[-IDs_zero,]#delete the gene has all 0.000 depression

for(i in 1:length(expression3)) {
  expression3[is.na(expression3[, i]), i] <- mean(expression3[, i], na.rm = T)
}#replace NA to mean of its sample expression


#remove batch effect---------------------------------------------------------------------------------
batch1=rep(1,length(Project1_id))
batch2=rep(2,length(Project2_id))
batch=c(batch1,batch2)
expression3=as.matrix(expression3)
combat_edata = ComBat(dat=expression3, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
write.table(combat_edata,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_CTLA4_removed_batch_expression.txt",quote = FALSE,row.names = TRUE,col.names = TRUE)

#make difference----------------------------------------------------------------------------------------
dplyr::select(as.data.frame(combat_edata),non_response$Run) ->non_response_expression
dplyr::select(as.data.frame(combat_edata),response$Run) ->response_expression
all_expression=cbind(response_expression,non_response_expression)
#orderthe expression profile by response and nonresponse to t.test

avg.R=apply(all_expression,1,function(x) median(x[1:nrow(response)]))
avg.NR=apply(all_expression,1,function(x) median(x[(nrow(response)+1):length(all_expression)]))
diff.avg=apply(all_expression,1,function(x) (median(x[1:nrow(response)])-median(x[(nrow(response)+1):length(all_expression)])))
p_value=apply(all_expression,1,function(x) t.test(x[1:nrow(response)],x[(nrow(response)+1):length(all_expression)])$p.value)
result=as.data.frame(cbind(all_expression,avg.R,avg.NR,diff.avg,p_value))

result=cbind(rownames(result),result)
rownames(result)=NULL
colnames(result)[1]="ensembl_ID"
write.table(result[,c(1,(length(result)-3),(length(result)-2),(length(result)-1),length(result))],"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/melanoma_CTLA4_DEG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#write all genes' difference

#filter significant different gene-------------------------------------
dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
  dplyr::filter(diff.avg>=2) -> up#57

dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
  dplyr::filter(diff.avg<=-2) -> down#105

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,up,by.x="EnsemblId",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(EnsemblId %in% up$ensembl_ID) ->up2
merge(relationship,down,by.x="EnsemblId",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(EnsemblId %in% down$ensembl_ID) ->down2
write.table(up2,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/up.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(down2,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/down.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#GO enrich for down-gene---------------------------------------------------------------------------------------
enrichGO(gene = up2$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,readable = TRUE)  ->ego_up
dotplot(ego_up)
enrichGO(gene = down2$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,readable = TRUE)->ego_down
dotplot(ego_down)
#write.table(as.data.frame(ego_up),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/up_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#0
write.table(as.data.frame(ego_down),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/down_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#1

#KEGG enrichment
enrichKEGG(gene=up2$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH") ->ekegg_up#0
enrichKEGG(gene=down2$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH") ->ekegg_down#2
browseKEGG(ekegg_down, 'hsa04110')

write.table(down_enrichKEGG,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/down_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
