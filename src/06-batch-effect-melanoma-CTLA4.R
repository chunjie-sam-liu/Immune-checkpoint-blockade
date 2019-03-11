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


#MDS test---------------------------------------------------------------------------------
#MDS before
a1=as.matrix(dist(t(expression3)))
voles.mds=cmdscale(a1,k=13,eig=T)
sum(abs(voles.mds$eig[1:2]))/sum(abs(voles.mds$eig))
sum((voles.mds$eig[1:2])^2)/sum((voles.mds$eig)^2)
x = voles.mds$points[,1]
y = voles.mds$points[,2]
point_color=c(rep("P1",length(Project1_id)),rep("P2",length(Project2_id)))
p=ggplot(data.frame(x,y,point_color),aes(x,y,label = colnames(a1),colour=point_color))+
  geom_point(shape=16,size=1)
ggsave(
  filename = 'MDS_before.pdf',
  plot = p,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/batch_effect/CTLA4',
  width = 6,
  height = 6.8
)

#MDS after removing batch effect with ComBat
a2=as.matrix(dist(t(combat_edata)))
voles.mds_2=cmdscale(a2,k=13,eig=T)
sum(abs(voles.mds_2$eig[1:2]))/sum(abs(voles.mds_2$eig))
sum((voles.mds_2$eig[1:2])^2)/sum((voles.mds_2$eig)^2)
x = voles.mds_2$points[,1]
y = voles.mds_2$points[,2]
point_color=c(rep("P1",length(Project1_id)),rep("P2",length(Project2_id)))
p=ggplot(data.frame(x,y,point_color),aes(x,y,label = colnames(a2),colour=point_color))+
  geom_point(shape=16,size=1)
ggsave(
  filename = 'MDS_ComBat.pdf',
  plot = p,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/batch_effect/CTLA4',
  width = 6,
  height = 6.8
)

#PCA test----------------------------------------------------------------------------------
#before
pca_result <- prcomp(t(expression3),scale=T) 
point_color=c(rep("P1",length(Project1_id)),rep("P2",length(Project2_id)))
ggplot(data.frame(pca_result$x[,1:2],point_color),aes(pca_result$x[,1],pca_result$x[,2],label = colnames(expression3),colour=point_color)) + geom_point(shape=16,size=1)->PCA_plot
ggsave(
  filename = 'PCA_before.pdf',
  plot = PCA_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/batch_effect/CTLA4',
  width = 6,
  height = 6.8
)

#combat
pca_result <- prcomp(t(combat_edata),scale=T) 
ggplot(data.frame(pca_result$x[,1:2],point_color),aes(pca_result$x[,1],pca_result$x[,2],label = colnames(combat_edata),colour=point_color)) + geom_point(shape=16,size=1)->PCA_plot
ggsave(
  filename = 'PCA_ComBat.pdf',
  plot = PCA_plot,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/batch_effect/CTLA4',
  width = 6,
  height = 6.8
)
