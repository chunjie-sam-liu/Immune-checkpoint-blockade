library(GSVA)
library(readxl)
library(magrittr)
library(GSEABase)

#load gene set,metadata and gene symbol,ensenml,entrezID relationship---------------------------------------
file_path = "/data/liucj/project/06-autophagy/GSEA/GSEA-collections-gmt"
genesets <- getGmt(paste(file_path,"C2_CURATED.gmt",sep="/"))

C6_ONCOGENIC_SIGNATURES.gmt
C7_IMMUNOLOGIC_SIGNATURES.gmt

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") ->metadata 
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE,sep="\t") -> relationship

#repeat SRP070710's result--------------------------------------------------------------------------------

read.table("/data/liull/immune-checkpoint-blockade/expression/all_FPKM_expression.txt",sep="\t",header = T,as.is = TRUE) %>% 
  dplyr::filter(gene_id %in% relationship$EnsemblId) %>% 
  merge(relationship[,2:3],.,by.x="EnsemblId",by.y="gene_id") %>%
  dplyr::select(-EnsemblId)->all_expression  #load all fo the original FPKM profile
 
metadata %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(SRA_Study=="SRP070710") %>%
  dplyr::select(Run,Response) -> SRP070710
dplyr::select(all_expression,c("Symbol",as.character(as.matrix(SRP070710)[,1]))) -> SRP070710_expression
rownames(SRP070710_expression) <- SRP070710_expression[,1]
SRP070710_expression <- SRP070710_expression[,-1]  #filter SRP070710 profile
  
SRP070710_GSVA_score=gsva(as.matrix(SRP070710_expression), genesets, min.sz=1, max.sz=999999, method="zscore",kcdf="Gaussian", abs.ranking=FALSE, verbose=TRUE)  # make GSVA

dplyr::filter(SRP070710,Response %in% c("CR","PR","R"))$Run %>%
  dplyr::select(as.data.frame(SRP070710_GSVA_score),.) -> SRP070710_response_gsva
dplyr::filter(SRP070710,Response %in% c("SD","PD","NR"))$Run %>%
  dplyr::select(as.data.frame(SRP070710_GSVA_score),.) -> SRP070710_non_response_gsva
SRP070710_ordered_GSVA=cbind(SRP070710_response_gsva,SRP070710_non_response_gsva) #reorder the SRP070710 GSVA score matrix

m=ncol(SRP070710_response_gsva)
n=ncol(SRP070710_ordered_GSVA)
SRP070710_avg_R=apply(SRP070710_ordered_GSVA,1,function(x) median(x[1:m]))
SRP070710_avg_NR=apply(SRP070710_ordered_GSVA,1,function(x) median(x[(m+1):n]))
SRP070710_diffAvg=apply(SRP070710_ordered_GSVA,1,function(x) abs((median(x[1:m])-median(x[(m+1):n]))/median(x[(m+1):n])))
SRP070710_p_value=apply(SRP070710_ordered_GSVA,1,function(x) t.test(x[1:m],x[(m+1):n])$p.value)
SRP070710_adj_p=p.adjust(SRP070710_p_value, "fdr")
SRP070710_result=as.data.frame(cbind(SRP070710_ordered_GSVA,SRP070710_avg_R,SRP070710_avg_NR,SRP070710_diffAvg,SRP070710_p_value,SRP070710_adj_p))

cbind(rownames(SRP070710_result),SRP070710_result) %>%
  dplyr::filter(SRP070710_diffAvg>=0.1) %>%
  dplyr::filter(SRP070710_p_value<=0.1)-> SRP070710_sig_sets
heatmap(as.matrix(SRP070710_sig_sets[,-c(1,30,31,32,33,34)]),hclustfun = hclust,Colv=NA,col=topo.colors(100))




#load removed batch effect's expression profile--------------------------------------------------------

read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_PD1_removed_batch_expression.txt",header = T,as.is = TRUE) %>% 
  cbind(rownames(.),.)->all_expression
colnames(all_expression)[1]="ensembl_ID"
all_expression %>%
  dplyr::filter(ensembl_ID %in% relationship$EnsemblId) %>%
  merge(relationship[,2:3],.,by.x="EnsemblId",by.y="ensembl_ID") %>%
  dplyr::select(-EnsemblId)->expression
rownames(expression) <- expression[,1]
expression <- expression[,-1] 

#order response or not
metadata %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(Run,Response) ->melanoma_PD1


#make GSVA score for all
GSVA_score <- gsva(data.matrix(expression), genesets, min.sz=1, max.sz=999999, method="zscore",kcdf="Gaussian", abs.ranking=FALSE, verbose=TRUE)

dplyr::filter(melanoma_PD1,Response %in% c("CR","PR","R","PRCR"))$Run ->response_ids
dplyr::filter(melanoma_PD1,Response %in% c("SD","PD","NR"))$Run ->non_response_ids
dplyr::select(as.data.frame(GSVA_score),response_ids,non_response_ids) ->ordered_GSVA

#IDs=as.vector(NULL)
#for (i in 1:nrow(ordered_GSVA)) {
#  x=as.character(ordered_GSVA[i,])
#  if (var(x)==0){IDs=c(IDs,i)}
#}
#ordered_GSVA=ordered_GSVA[-IDs,]#delete the gene sets hava too many NA expression's gene

avg.R=apply(ordered_GSVA,1,function(x) median(x[1:length(response_ids)]))
avg.NR=apply(ordered_GSVA,1,function(x) median(x[(length(response_ids)+1):ncol(ordered_GSVA)]))
diff.avg=apply(ordered_GSVA,1,function(x) (median(x[1:length(response_ids)])-median(x[(length(response_ids)+1):ncol(ordered_GSVA)])))
p_value=apply(ordered_GSVA,1,function(x) t.test(x[1:length(response_ids)],x[(length(response_ids)+1):ncol(ordered_GSVA)])$p.value)
result=as.data.frame(cbind(ordered_GSVA,avg.R,avg.NR,diff.avg,p_value))

cbind(rownames(result),result) %>%
  dplyr::filter(abs(diff.avg)>=0.1) %>%
  dplyr::filter(p_value<=0.01)-> sig_sets
write.table(sig_sets,"/data/liull/immune-checkpoint-blockade/different_expression/melanoma/sig_sets.txt",sep="\t",quote=FALSE,col.names = TRUE,row.names = FALSE)

rownames(sig_sets) = sig_sets[,1]
dplyr::select(sig_sets,2:(ncol(sig_sets)-4))->a


heatmap(as.matrix(a), Colv=NA,ColSideColors=c(rep("purple", 44), rep("orange", 117)),col=colorRampPalette(c("green", "black","red"))(256),cexRow = 0.5,cexCol = 0.2)
