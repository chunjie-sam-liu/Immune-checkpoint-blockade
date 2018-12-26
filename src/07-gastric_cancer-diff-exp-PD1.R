library(magrittr)
library(readr)
library(readxl)
library(dplyr)
library(sva)

#筛选出gastric cancer的RNA-seq的anti-PD1的数据的应答信息-------------
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="gastric cancer") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

#select response and non_response's sample id and project id-------------------
rbind(dplyr::filter(metadata,Response=="CR"),dplyr::filter(metadata,Response=="PR")) ->response
rbind(dplyr::filter(metadata,Response=="SD"),dplyr::filter(metadata,Response=="PD")) ->non_response

#make difference----------------------------
read.table("/data/liull/immune-checkpoint-blockade/expression/all_expression.txt",sep="\t",header = T,as.is = TRUE) ->data
dplyr::select(data,gene_id,response$Run) ->response_expression
dplyr::select(data,non_response$Run) ->non_response_expression
all_expression=cbind(response_expression,non_response_expression)
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
response_Mean=apply(all_expression3,1,function(x) mean(x[1:nrow(response)]))
non_response_Mean=apply(all_expression3,1,function(x) mean(x[(nrow(response)+1):length(all_expression3)]))
FC=apply(all_expression3,1,function(x) (mean(x[1:nrow(response)])+0.01)/(mean(x[(nrow(response)+1):length(all_expression3)])+0.01))
p_value=apply(all_expression3,1,function(x) t.test(x[1:nrow(response)],x[(nrow(response)+1):length(all_expression3)])$p.value)
result=as.data.frame(cbind(all_expression3,response_Mean,non_response_Mean,FC,p_value))

result=cbind(rownames(result),result)
rownames(result)=NULL
colnames(result)[1]="ensembl_ID"
result$ensembl_ID <- as.character(result$ensembl_ID)

dplyr::filter(as.data.frame(result),p_value<=0.1) %>%
  dplyr::filter(FC>=2) -> up
dplyr::filter(as.data.frame(result),p_value<=0.1) %>%
  dplyr::filter(FC<=0.5) -> down

write.table(result,"/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/all_values.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(up,"/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/up_values.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(down,"/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/down_values.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
