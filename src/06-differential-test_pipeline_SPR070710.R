library(readxl)
library(dplyr)
library(magrittr)
library(NOISeq)
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(SRA_Study=="SRP070710") %>%
  #dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(Run,Response) ->metadata

dplyr::filter(metadata,Response %in% c("CR","PR")) ->response
dplyr::filter(metadata,Response=="PD") ->non_response

#select response and non_response
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
read.table("/data/liull/immune-checkpoint-blockade/expression/all_FPKM_expression_2.txt",sep="\t",header = T,as.is = TRUE) ->data1

data1%>%
  dplyr::select(gene_id,response$Run,non_response$Run)%>%
  dplyr::filter(gene_id %in% grep("ENSG",data1$gene_id,value=T))->SRP070710_expr
# %>%
#   dplyr::filter(gene_id %in% relationship$Ensembl_ID) %>%
#   merge(relationship,.,by.x="Ensembl_ID",by.y="gene_id")%>%
#   dplyr::select(-Ensembl_ID,-GeneID)->SRP070710_expr  #[1] 21003    29

# factors=factor(SRP070710_expr$Symbol)
# merged_expression=tapply(SRP070710_expr[,2],factors,median)
# for (i in 3:ncol(SRP070710_expr)) {
#   temp=tapply(SRP070710_expr[,i],factors,median)
#   merged_expression=cbind(merged_expression,temp)
# }
# colnames(merged_expression)=colnames(SRP070710_expr)[2:ncol(SRP070710_expr)]  #trans ensembl id to symbol and merged
# dim(merged_expression)#[1] 20986    28

# keep <- rowSums(merged_expression>3,na.rm = TRUE) > 0
# merged_expression[keep,]->filtered_expr  #[1] 11165    28
rownames(SRP070710_expr)=SRP070710_expr$gene_id
SRP070710_expr=SRP070710_expr[,-1]
SRP070710_expr[is.na(SRP070710_expr)]<-0

median_R=apply(SRP070710_expr, 1, function(x) median(x[1:15]))
median_NR=apply(SRP070710_expr, 1, function(x) median(x[16:28]))
diff_median=apply(SRP070710_expr, 1, function(x) (median(x[1:15])+0.001)/(median(x[16:28])+0.001))
p.val=apply(SRP070710_expr, 1, function(x) wilcox.test(x[1:15],x[16:28],alternative="two.sided")$p.value)

cbind(median_R,median_NR,diff_median,p.val)%>%
  as.data.frame(stringsAsFactors=FALSE)->result
cbind(rownames(SRP070710_expr),result,stringsAsFactors=FALSE)->result

dplyr::filter(result,p.val<0.1)%>%
  dplyr::filter(diff_median>2)->up
dplyr::filter(result,p.val<0.1)%>%
  dplyr::filter(diff_median< 0.5)->down


readxl::read_excel("/data/liull/immune-checkpoint-blockade/different_expression/test_SRP070710/1-s2.0-S009286741630215X-mmc2.xlsx",col_names = TRUE,sheet="S2A",skip = 2) ->standard_diff_genes
dplyr::filter(standard_diff_genes,diffAvg>0)%>%
  merge(relationship,.,by.x="Symbol",by.y="Gene")%>%
  dplyr::select(Ensembl_ID)->standard_up

standard_down=dplyr::filter(standard_diff_genes,diffAvg<0)%>%
  merge(relationship,.,by.x="Symbol",by.y="Gene")%>%
  dplyr::select(Ensembl_ID)->standard_down

length(intersect(standard_up$Ensembl_ID,up$`rownames(SRP070710_expr)`))
# median(R)-median(NR): 30  FC:59
length(intersect(standard_down$Ensembl_ID,down$`rownames(SRP070710_expr)`))
# median(R)-median(NR): 94  FC:179
 
venn.diagram(list(SRP094781=up$`rownames(SRP070710_expr)`, SRP094781_standard=standard_up$Ensembl_ID),
             cex=2,margin = 0.2,imagetype = "png",
             filename="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP094781_one_method_up.png",
             fill=c("red","yellow"),na="remove")
venn.diagram(list(SRP094781=down$`rownames(SRP070710_expr)`, SRP094781_standard=standard_down$Ensembl_ID),
             cex=2,margin = 0.2,imagetype = "png",
             filename="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP094781_one_method_down.png",
             fill=c("red","yellow"),na="remove")
