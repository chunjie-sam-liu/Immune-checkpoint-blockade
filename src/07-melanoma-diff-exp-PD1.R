library(magrittr)
library(readr)
library(readxl)
library(dplyr)
library(sva)

#筛选出黑色素瘤的RNA-seq的anti-PD1的数据的应答信息
readxl::read_excel("/data/liull/immune-checkpoint-blockade/04-all-metadata.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

#select response and non_response's sample id and project id
CR=dplyr::filter(metadata,Response=="CR")
PR=dplyr::filter(metadata,Response=="PR")
R=dplyr::filter(metadata,Response=="R")
rbind(CR,PR) %>%
  rbind(R) ->response
SD=dplyr::filter(metadata,Response=="SD")
PD=dplyr::filter(metadata,Response=="PD")
NR=dplyr::filter(metadata,Response=="NR")
rbind(SD,PD) %>%
  rbind(NR) ->non_response

#expression prepare for batch effect
read.table("/data/liull/immune-checkpoint-blockade/expression/all_expression.txt",sep="\t",header = T,as.is = TRUE) ->data1
Project=unique(metadata$SRA_Study)
dplyr::filter(metadata,SRA_Study==Project[1]) %>%
  dplyr::select(Run)  %>%
  as.matrix() %>%
  as.character()->Project1_id
dplyr::select(data1,gene_id,Project1_id) ->Project1_expression

dplyr::filter(metadata,SRA_Study==Project[2]) %>%
  dplyr::select(Run)  %>%
  as.matrix() %>%
  as.character()->Project2_id
dplyr::select(data1,Project2_id) ->Project2_expression
expression1=cbind(Project1_expression,Project2_expression)
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


#remove batch effect
batch1=rep(1,length(Project1_id))
batch2=rep(2,length(Project2_id))
batch=c(batch1,batch2)
expression3=as.matrix(expression3)
combat_edata = ComBat(dat=expression3, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

#make difference
dplyr::select(as.data.frame(combat_edata),non_response$Run) ->non_response_expression
dplyr::select(as.data.frame(combat_edata),response$Run) ->response_expression
all_expression=cbind(response_expression,non_response_expression)
#orderthe expression profile by response and nonresponse to t.test

response_Mean=apply(all_expression,1,function(x) mean(x[1:nrow(response)]))
non_response_Mean=apply(all_expression,1,function(x) mean(x[(nrow(response)+1):length(all_expression)]))
FC=apply(all_expression,1,function(x) (mean(x[1:nrow(response)])+0.01)/(mean(x[(nrow(response)+1):length(all_expression)])+0.01))
p_value=apply(all_expression,1,function(x) t.test(x[1:nrow(response)],x[(nrow(response)+1):length(all_expression)])$p.value)
result=as.data.frame(cbind(response_Mean,non_response_Mean,FC,p_value))

result=cbind(rownames(result),result)
rownames(result)=NULL
colnames(result)[1]="ensembl_ID"
result$ensembl_ID <- as.character(result$ensembl_ID)

dplyr::filter(as.data.frame(result),p_value<=0.1) %>%
  dplyr::filter(FC>=2) -> up
dplyr::filter(as.data.frame(result),p_value<=0.1) %>%
  dplyr::filter(FC<=0.5) -> down
write.table(result,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/all_values.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(up,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/up_values.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(down,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/down_values.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
