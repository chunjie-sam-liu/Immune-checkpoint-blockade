library(clusterProfiler)
library(org.Hs.eg.db)
library(clusterProfiler)
library(readxl)
library(magrittr)


#GSVA
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE,sep="\t") -> relationship

#melanoma_PD1----------------------------------------------------------------------------------------

read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_PD1_removed_batch_expression.txt",header = T,as.is = TRUE) %>% 
  cbind(rownames(.),.)->all_expression
colnames(all_expression)[1]="ensembl_ID"


all_expression %>%
  dplyr::filter(ensembl_ID %in% relationship$EnsemblId) %>%
  merge(relationship[,2:3],.,by.x="EnsemblId",by.y="ensembl_ID") %>%
  dplyr::select(-EnsemblId)->expression
rownames(expression) <- expression[,1]
expression <- expression[,-1] 

file_path = "/data/liull/reference/GSEA-gmt"
gmt_2<- getGmt(paste(file_path,"c2.all.v6.2.symbols.gmt",sep="/"))
gmt_7<- getGmt(paste(file_path,"c7.all.v6.2.symbols.gmt",sep="/"))

GSVA_score_gmt2 <- gsva(data.matrix(expression), gmt_2, min.sz=1, max.sz=999999, method="zscore",kcdf="Gaussian", abs.ranking=FALSE, verbose=TRUE)
GSVA_score_gmt7 <- gsva(data.matrix(expression), gmt_7, min.sz=1, max.sz=999999, method="zscore",kcdf="Gaussian", abs.ranking=FALSE, verbose=TRUE)

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") ->metadata 
metadata %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(Run,Response) ->melanoma_PD1
dplyr::filter(melanoma_PD1,Response %in% c("CR","PR","R","PRCR"))$Run ->response_ids
dplyr::filter(melanoma_PD1,Response %in% c("SD","PD","NR"))$Run ->non_response_ids
dplyr::select(as.data.frame(GSVA_score_gmt2),response_ids,non_response_ids) ->ordered_GSVA

avg.R=apply(ordered_GSVA,1,function(x) median(x[1:length(response_ids)]))
avg.NR=apply(ordered_GSVA,1,function(x) median(x[(length(response_ids)+1):ncol(ordered_GSVA)]))
diff.avg=apply(ordered_GSVA,1,function(x) (median(x[1:length(response_ids)])-median(x[(length(response_ids)+1):ncol(ordered_GSVA)])))
p_value=apply(ordered_GSVA,1,function(x) t.test(x[1:length(response_ids)],x[(length(response_ids)+1):ncol(ordered_GSVA)])$p.value)
diff_GSVA_gmt2=as.data.frame(cbind(ordered_GSVA,avg.R,avg.NR,diff.avg,p_value))

cbind(rownames(diff_GSVA_gmt2),diff_GSVA_gmt2) %>%
  dplyr::filter(abs(diff.avg)>=1) %>%
  dplyr::filter(p_value<=0.05)-> sig_sets
write.table(sig_sets,"/data/liull/immune-checkpoint-blockade/GSVA/melanoma_PD1/sig_sets.txt",sep="\t",quote=FALSE,col.names = TRUE,row.names = FALSE)
