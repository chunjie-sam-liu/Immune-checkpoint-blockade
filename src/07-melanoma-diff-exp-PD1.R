library(magrittr)
library(readr)
library(readxl)
library(dplyr)
library(sva)

#筛选出黑色素瘤的RNA-seq的anti-PD1的数据的应答信息-------------------------------------------------------------
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

#select response and non_response's sample id and project id
CR=dplyr::filter(metadata,Response %in% c("CR","PR","PRCR","R")) -> response
SD=dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response

#expression prepare for batch effect----------------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/expression/all_FPKM_expression.txt",sep="\t",header = T,as.is = TRUE) ->data1
Project=unique(metadata$SRA_Study)
dplyr::filter(metadata,SRA_Study==Project[1]) %>%
  dplyr::select(Run)  %>%
  as.matrix() %>%
  as.character()->Project1_id

dplyr::filter(metadata,SRA_Study==Project[2]) %>%
  dplyr::select(Run)  %>%
  as.matrix() %>%
  as.character()->Project2_id

dplyr::filter(metadata,SRA_Study==Project[3]) %>%
  dplyr::select(Run)  %>%
  as.matrix() %>%
  as.character()->Project3_id

expression1=dplyr::select(data1,gene_id,Project1_id,Project2_id,Project3_id)
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
batch3=rep(3,length(Project3_id))
batch=c(batch1,batch2,batch3)
expression3=as.matrix(expression3)
combat_edata = ComBat(dat=expression3, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
write.table(combat_edata,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_PD1_removed_batch_expression.txt",quote = FALSE,row.names = TRUE,col.names = TRUE)

#make difference----------------------------------------------------------------------------------------
dplyr::select(as.data.frame(combat_edata),non_response$Run) ->non_response_expression
dplyr::select(as.data.frame(combat_edata),response$Run) ->response_expression
all_expression=cbind(response_expression,non_response_expression)
#orderthe expression profile by response and nonresponse to t.test

avg.R=apply(all_expression,1,function(x) median(x[1:nrow(response)]))
avg.NR=apply(all_expression,1,function(x) median(x[(nrow(response)+1):length(all_expression)]))
diff.avg=apply(all_expression,1,function(x) (median(x[1:nrow(response)])-median(x[(nrow(response)+1):length(all_expression)])))
p_value=apply(all_expression,1,function(x) wilcox.test(x[1:nrow(response)],x[(nrow(response)+1):length(all_expression)])$p.value)
result=as.data.frame(cbind(all_expression,avg.R,avg.NR,diff.avg,p_value))

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
result=cbind(rownames(result),result)
colnames(result)[1]="ensembl_ID"
merge(relationship,result,by.x="EnsemblId",by.y="ensembl_ID")->result
write.table(result,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_PD1_DEG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#filter significant different gene(same cutoff with the-28-sample-tableS2:0.1 , 1)-------------------------------------


dplyr::filter(as.data.frame(result),p_value<=0.1) %>%
  dplyr::filter(diff.avg>=1) -> up#232
rownames(up)=up[,2]

dplyr::filter(as.data.frame(result),p_value<=0.1) %>%
  dplyr::filter(diff.avg<=-1) -> down#699
rownames(down)=down[,3]


#heatmap for down-gene(p<0.01)--------------------------------------------------------------------------------
dplyr::filter(down,p_value<=0.01) %>%
  dplyr::filter(diff.avg<=-2)-> down_for_map
rownames(down_for_map)=down_for_map[,3]
heatmap(as.matrix(down_for_map[,4:(ncol(down_for_map)-4)]),Colv=NA,ColSideColors=c(rep("purple", 44), rep("orange", 117)),col=colorRampPalette(c("green", "black","red"))(256),cexRow = 0.3,cexCol = 0.2)

write.table(down_for_map,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/down_for_heatmap.txt",sep="\t",quote=FALSE,col.names = TRUE,row.names = FALSE)
#GO enrich for down-gene---------------------------------------------------------------------------------------
enrichGO(gene = down_for_map$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE) %>%
  as.data.frame()->down_for_map_enrich
