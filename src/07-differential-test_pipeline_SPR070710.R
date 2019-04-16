library(readxl)
library(dplyr)
library(magrittr)
library(NOISeq)
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(SRA_Study=="SRP070710") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(Run,Response) ->metadata
#filter SRP070710  metadata
dplyr::filter(metadata,Response %in% c("CR","PR")) ->response
dplyr::filter(metadata,Response=="PD") ->non_response

#select response and non_response

read.table("/data/liull/immune-checkpoint-blockade/expression/all_FPKM_expression_2.txt",sep="\t",header = T,as.is = TRUE) ->data1

dplyr::select(data1,gene_id,response$Run,non_response$Run)->all_expression
rownames(all_expression)=all_expression[,1]
all_expression=all_expression[,-1]


all_expression[is.na(all_expression)]<-0
# Sum_zero=apply(all_expression,1,function(x) sum(as.numeric(x),na.rm = TRUE))
# IDs_zero=which(Sum_zero==0)
# all_expression_2=all_expression[-IDs_zero,]#delete the gene has all 0.000 expression


b=nrow(response)
c=nrow(response)+1
d=ncol(all_expression)
Pvalue<-c(rep(0,nrow(all_expression)))
FC<-c(rep(0,nrow(all_expression)))
for(i in 1:nrow(all_expression)){
  
  if(sd(all_expression[i,1:b],na.rm = TRUE)==0){
    Pvalue[i] <- "NA"
    FC[i]<- "NA"
  }else{
    y=wilcox.test(as.numeric(all_expression[i,1:b]),as.numeric(all_expression[i,c:d]),alternative = "two.sided")
    Pvalue[i]<-y$p.value
    FC[i]<-median(as.numeric(all_expression[i,1:b]))-median(as.numeric(all_expression[i,c:d]))
  }

}

avg_response=apply(all_expression[,1:b],1,median)
avg_nonresponse=apply(all_expression[,c:d],1,median)

out<-cbind(all_expression,avg_response,avg_nonresponse,Pvalue,FC)
out$Pvalue=as.numeric(as.character(out$Pvalue))
out$FC=as.numeric(as.character(out$FC))
# > dim(out)
# [1] 49325    31

del=which(is.na(out$Pvalue)==TRUE)
out2=out[-del,]#delete the Pvalue is NA 
# > dim(out2)
# [1] 44128    31

out2=cbind(rownames(out2),out2)
rownames(out2)=NULL
colnames(out2)[1]="gene_id"
out2$gene_id <- as.character(out2$gene_id)#change rownames to out[,1]
#write.table(out2,"/data/liull/immune-checkpoint-blockade/different_expression/test_SRP070710/all_diff_expression.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

dplyr::filter(out2,Pvalue<=0.1) -> result
dplyr::filter(result,FC>=1) -> up_in_response
dplyr::filter(result,FC<= -1) -> down_in_response
# write.table(up_in_response,"/data/liull/immune-checkpoint-blockade/different_expression/test_SRP070710/up_0.1_2.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
# write.table(down_in_response,"/data/liull/immune-checkpoint-blockade/different_expression/test_SRP070710/down_0.1_2.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)



readxl::read_excel("/data/liull/immune-checkpoint-blockade/different_expression/test_SRP070710/1-s2.0-S009286741630215X-mmc2.xlsx",col_names = TRUE,sheet="S2A",skip = 2) ->standard_diff_genes
standard_up=dplyr::filter(standard_diff_genes,diffAvg>0)
standard_down=dplyr::filter(standard_diff_genes,diffAvg<0)

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI_2.txt",sep="\t",header = T,as.is = TRUE) ->relationship

dplyr::filter(up_in_response,gene_id %in% relationship$EnsemblId) %>%
  merge(relationship,.,by.x="EnsemblId",by.y="gene_id")%>%
  dplyr::select(Symbol)->up_in_response_symbol
length(dplyr::intersect(standard_up$Gene,up_in_response_symbol$Symbol))#30

dplyr::filter(down_in_response,gene_id %in% relationship$EnsemblId) %>%
  merge(relationship,.,by.x="EnsemblId",by.y="gene_id")%>%
  dplyr::select(Symbol)->down_in_response_symbol
length(dplyr::intersect(standard_down$Gene,down_in_response_symbol$Symbol))#76



# ID_up=as.character(as.data.frame(standard_up[,1])$Gene)
# ID_down=as.character(as.data.frame(standard_down[,1])$Gene)
# read.table("/data/liull/salmon_tutorial/Homo_quants/gene_transcript_relationship.txt",header = T) -> relationship
# list=relationship[,c(1,3)]
# ID_list_up=list[match(ID_up,list[,"gene_name"]),]#从gene_transcript_relationship.txt中提取出up的基因的ensembl id。
# ID_list_down=list[match(ID_down,list[,"gene_name"]),]#从gene_transcript_relationship.txt中提取出down的基因的ensembl id。
# 
# length(dplyr::intersect(ID_list_up$gene_id,up_in_response$gene_id))
# length(dplyr::intersect(ID_list_down$gene_id,down_in_response$gene_id))


