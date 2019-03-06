library(GSVA)
library(readxl)
library(magrittr)
library(GSEABase)

#load gene set and gene symbol,ensenml,entrezID relationship---------------------------------------
file_path = "/data/liull/reference/GSEA-gmt"
genesets_c7<- getGmt(paste(file_path,"c7.all.v6.2.symbols.gmt",sep="/"))

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE,sep="\t") -> relationship

#load melanoma_PD1_removed_batch_expression and translate Ensembl id to symbol--------------------------------------------

read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_PD1_removed_batch_expression.txt",header = T,as.is = TRUE) %>% 
  cbind(rownames(.),.)->all_expression
colnames(all_expression)[1]="ensembl_ID"
all_expression %>%
  dplyr::filter(ensembl_ID %in% relationship$EnsemblId) %>%
  merge(relationship[,2:3],.,by.x="EnsemblId",by.y="ensembl_ID") %>%
  dplyr::select(-EnsemblId)->expression
rownames(expression) <- expression[,1]
expression <- expression[,-1] 

#make GSVA score for melanoma_PD1
GSVA_score_c7 <- gsva(data.matrix(expression), genesets_c7, min.sz=1, max.sz=999999, method="zscore",kcdf="Gaussian", abs.ranking=FALSE, verbose=TRUE)

#order response or not
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(Run,Response) ->melanoma_PD1
dplyr::filter(melanoma_PD1,Response %in% c("CR","PR","R","PRCR"))$Run ->response_ids
dplyr::filter(melanoma_PD1,Response %in% c("SD","PD","NR"))$Run ->non_response_ids

dplyr::select(as.data.frame(GSVA_score_c7),response_ids,non_response_ids) ->ordered_GSVAscore_c7

#IDs=as.vector(NULL)
#for (i in 1:nrow(ordered_GSVA)) {
#  x=as.character(ordered_GSVA[i,])
#  if (var(x)==0){IDs=c(IDs,i)}
#}
#ordered_GSVA=ordered_GSVA[-IDs,]#delete the gene sets hava too many NA expression's gene

avg.R=apply(ordered_GSVAscore_c7,1,function(x) median(x[1:length(response_ids)]))
avg.NR=apply(ordered_GSVAscore_c7,1,function(x) median(x[(length(response_ids)+1):ncol(ordered_GSVAscore_c7)]))
diff.avg=apply(ordered_GSVAscore_c7,1,function(x) (median(x[1:length(response_ids)])-median(x[(length(response_ids)+1):ncol(ordered_GSVAscore_c7)])))
p_value=apply(ordered_GSVAscore_c7,1,function(x) t.test(x[1:length(response_ids)],x[(length(response_ids)+1):ncol(ordered_GSVAscore_c7)])$p.value)
result=as.data.frame(cbind(ordered_GSVAscore_c7,avg.R,avg.NR,diff.avg,p_value))

cbind(rownames(result),result) %>%
  dplyr::filter(abs(diff.avg)>=0.1) %>%
  dplyr::filter(p_value<=0.05)-> sig_sets
write.table(sig_sets,"/data/liull/immune-checkpoint-blockade/GSVA/melanoma_PD1/c7_sig_sets.txt",sep="\t",quote=FALSE,col.names = TRUE,row.names = FALSE)

#heatmap
rownames(sig_sets) = sig_sets[,1]
dplyr::select(sig_sets,2:(ncol(sig_sets)-4))->a

annotation_col = data.frame(SampleClass = factor(rep(c("response", "non-response"), c(length(response_ids),length(non_response_ids)))))
rownames(annotation_col)=colnames(a)
pheatmap(a, annotation_col = annotation_col,cluster_cols = FALSE,scale="column")

heatmap(as.matrix(a), Colv=NA,ColSideColors=c(rep("purple", 44), rep("orange", 117)),col=colorRampPalette(c("green", "black","red"))(256),cexRow = 0.5,cexCol = 0.3)->p
