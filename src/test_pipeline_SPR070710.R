library(readxl)
library(dplyr)
library(magrittr)
library(NOISeq)
readxl::read_excel("/data/liull/immune-checkpoint-blockade/04-all-metadata.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(SRA_Study=="SRP070710") %>%
  dplyr::select(Run,Response) ->metadata
#filter SRP070710  metadata
CR=dplyr::filter(metadata,Response=="CR")
PR=dplyr::filter(metadata,Response=="PR")
SD=dplyr::filter(metadata,Response=="SD")
rbind(CR,PR) %>%
  rbind(SD) ->response

dplyr::filter(metadata,Response=="PD") ->non_response

#select response and non_response

expression_files = read.table("/data/liull/immune-checkpoint-blockade/expression/all_expression.txt",sep="\t",header = T,as.is = TRUE) ->data1
rownames(data1)=data1[,1]
data1=data1[,-1]

dplyr::select(data1,one_of(response$Run)) ->response_expression
dplyr::select(data1,one_of(non_response$Run)) ->non_response_expression
all_expression=cbind(response_expression,non_response_expression)
all_expression[is.na(all_expression)]<-0

b=ncol(response_expression)
c=ncol(response_expression)+1
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
    FC[i]<-(mean(as.numeric(all_expression[i,1:b]))+0.001)/(mean(as.numeric(all_expression[i,c:d]))+0.001)
  }

}
fdr=p.adjust(Pvalue, "fdr")
out<-cbind(all_expression,Pvalue,FC,fdr)
out$Pvalue=as.numeric(as.character(out$Pvalue))
out$FC=as.numeric(as.character(out$FC))


del=which(is.na(out$Pvalue)==TRUE)
out2=out[-del,]#delete the Pvalue is NA 

out2=cbind(rownames(out2),out2)
rownames(out2)=NULL
colnames(out2)[1]="gene_id"
out2$gene_id <- as.character(out2$gene_id)#change rownames to out[,1]

dplyr::filter(out2,Pvalue<=0.1) -> result
dplyr::filter(result1,FC>=2) -> up_in_response
dplyr::filter(result1,FC<=0.5) -> down_in_response

readxl::read_excel("/data/liull/immune-checkpoint-blockade/different_expression/1-s2.0-S009286741630215X-mmc2.xlsx",col_names = TRUE,sheet="S2A",skip = 2) ->standard_diff_genes
standard_up=dplyr::filter(standard_diff_genes,diffAvg>0)
standard_down=dplyr::filter(standard_diff_genes,diffAvg<0)
ID_up=as.character(as.data.frame(standard_up[,1])$Gene)
ID_down=as.character(as.data.frame(standard_down[,1])$Gene)
read.table("/data/liull/salmon_tutorial/Homo_quants/gene_transcript_relationship.txt",col.names = TRUE) -> relationship
list=relationship[,c(1,3)]
ID_list_up=list[match(ID_up,list[,"gene_name"]),]
ID_list_down=list[match(ID_down,list[,"gene_name"]),]

length(dplyr::intersect(ID_list_up$gene_id,up_in_response$gene_id))
length(dplyr::intersect(ID_list_down$gene_id,down_in_response$gene_id))


