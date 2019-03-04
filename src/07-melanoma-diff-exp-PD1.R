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

a=1:nrow(expression3)
for(j in 1:length(expression3)){
  which(expression3[,j]<10) ->b
  intersect(a,b)->a
}
expression4=expression3[-a,]

for(i in 1:length(expression4)) {
  expression4[is.na(expression4[, i]), i] <- mean(expression4[, i], na.rm = T)
}#replace NA to mean of its sample expression




#remove batch effect---------------------------------------------------------------------------------
batch1=rep(1,length(Project1_id))
batch2=rep(2,length(Project2_id))
batch3=rep(3,length(Project3_id))
batch=c(batch1,batch2,batch3)
expression4=as.matrix(expression4)
combat_edata = ComBat(dat=expression4, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
write.table(combat_edata,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_PD1_removed_batch_expression.txt",quote = FALSE,row.names = TRUE,col.names = TRUE)

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
write.table(result,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/melanoma_PD1_DEG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#write all genes' difference


#filter significant different gene-------------------------------------
dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
  dplyr::filter(diff.avg>=2) -> up#54

dplyr::filter(as.data.frame(result),p_value<=0.05) %>%
  dplyr::filter(diff.avg<=-2) -> down#199

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,up,by.x="EnsemblId",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(EnsemblId %in% up$ensembl_ID) ->up2
merge(relationship,down,by.x="EnsemblId",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(EnsemblId %in% down$ensembl_ID) ->down2
write.table(up2,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/up.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(down2,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/down.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)


#GO enrichment--------------------------------------------------------------------------------------
enrichGO(gene = up2$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,readable = TRUE)->ego_up
dotplot(ego_up,showCategory=20)
enrichGO(gene = down2$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.05,readable = TRUE) ->ego_down
dotplot(ego_down,showCategory=20)

write.table(as.data.frame(ego_up),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/up_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#171
write.table(as.data.frame(ego_down),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/down_enrichGO.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#58

#KEGG enrichment
enrichKEGG(gene=up2$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH") ->ekegg_up#32
enrichKEGG(gene=down2$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH")->ekegg_down#1
browseKEGG(ekegg_up, 'hsa04934')
browseKEGG(ekegg_down, 'hsa04934')

write.table(as.data.frame(ekegg_up),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/up_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(as.data.frame(ekegg_down),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/down_enrichKEGG.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#Reactome enrichment------------------------------------------------------
enrichPathway(gene=up2$GeneID,pvalueCutoff=0.05, readable=T)->eReactome_up#17
dotplot(eReactome_up, showCategory=17)
enrichPathway(gene=down2$GeneID,pvalueCutoff=0.05, readable=T)->eReactome_down#11
dotplot(eReactome_down, showCategory=11)

write.table(as.data.frame(eReactome_up),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/up_enrichReactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(as.data.frame(eReactome_down),"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/down_enrichReactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)



#heatmap for down-gene(p<0.01)----------------------------------------------------------
#dplyr::filter(down,p_value<=0.01) %>%
#  dplyr::filter(diff.avg<=-2)-> down_for_map
#rownames(down_for_map)=down_for_map[,3]
#heatmap(as.matrix(down_for_map[,4:(ncol(down_for_map)-4)]),Colv=NA,ColSideColors=c(rep("purple", 44), rep("orange", 117)),col=colorRampPalette(c("green", "black","red"))(256),cexRow = 0.3,cexCol = 0.2)

#Heatmap(ComplexHeatmap)
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/up.txt",sep="\t",header = T,as.is = TRUE) ->up
rownames(up)=up[,1]
up[order(up$p_value,decreasing=T),]->up2#diff.avg ordered from large to small
up2[,4:(ncol(up2)-4)] -> up2
up2[1:10,] %>% apply(2,function(x) scale(x)) ->a
rownames(a)=rownames(up2)[1:10]
Heatmap(as.matrix(a),cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2))


c(rep("1",44),rep("0",117)) %>% as.matrix() ->label
up3=rbind(label[,1],up2)
up3=as.data.frame(t(up3))
annot_up <- data.frame(response_VS_none = up3$`1`)
colors = list(response_VS_none = c("1" = "green", "0" = "gray"))
ha <- HeatmapAnnotation(annot_up, col = colors)
Heatmap(as.matrix(a),column_names_gp = gpar(fontsize = 2),top_annotation = ha,cluster_columns = FALSE)#Heatmap for ten top P

rownames(down2)=down2$EnsemblId
down2[,-c(1,2,3,(ncol(down2)-3),(ncol(down2)-2),(ncol(down2)-1),(ncol(down2)))]->a
annotation_col = data.frame(GeneClass = factor(rep(c("response", "non-response"), c(length(response_expression),length(non_response_expression)))))
rownames(annotation_col)=colnames(a)
pheatmap(a, annotation_col = annotation_col,cluster_cols = FALSE,scale="column")

