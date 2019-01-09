library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(readxl)
library(magrittr)
library(GSEABase)

#load gene set------------------------------------------------------------------------------------------
file_path = "/data/liull/immune-checkpoint-blockade/enrichment/melanoma"
genesets <- getGmt(paste(file_path,"c2.cgp.v6.2.symbols.gmt",sep="/"))

#load expression profile---------------------------------------------------------------------------------
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE,sep="\t") -> relationship
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/all_values.txt",sep="\t",header = T,as.is = TRUE) %>% 
  dplyr::select(.,1:(ncol(.)-4))%>%
  dplyr::filter(ensembl_ID %in% relationship$EnsemblId) %>%
  merge(relationship[,2:3],.,by.x="EnsemblId",by.y="ensembl_ID") %>%
  dplyr::select(-EnsemblId)->expression
rownames(expression) <- expression[,1]
expression <- expression[,-1]   

#get GSVA score------------------------------------------------------------------------------------------
topMatrixGSVA <- gsva(data.matrix(expression), genesets, min.sz=1, max.sz=999999, kcdf="Poisson", abs.ranking=FALSE, verbose=TRUE)


#make difference-----------------------------------------------------------------------------------------
readxl::read_excel("/data/liull/immune-checkpoint-blockade/04-all-metadata.xlsx",col_names = TRUE,sheet="SRA") ->metadata 
metadata %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(Run,Response) ->melanoma_PD1

dplyr::filter(melanoma_PD1,Response %in% c("CR","PR","R"))$Run ->response_ids
dplyr::filter(melanoma_PD1,Response %in% c("SD","PD","NR"))$Run ->non_response_ids
dplyr::select(as.data.frame(topMatrixGSVA),response_ids) ->response_gsva
dplyr::select(as.data.frame(topMatrixGSVA),non_response_ids) ->non_response_gsva
ordered_GSVA=cbind(response_gsva,non_response_gsva)

IDs=as.vector(NULL)
for (i in 1:nrow(ordered_GSVA)) {
  x=as.character(ordered_GSVA[i,])
  if (var(x)==0){IDs=c(IDs,i)}
}
ordered_GSVA=ordered_GSVA[-IDs,]#delete the gene sets hava too many NA expression's gene

response_Mean=apply(ordered_GSVA,1,function(x) mean(x[1:length(response_ids)]))
non_response_Mean=apply(ordered_GSVA,1,function(x) mean(x[(length(response_ids)+1):nrow(melanoma_PD1)]))
FC=apply(ordered_GSVA,1,function(x) (mean(x[1:length(response_ids)])+0.01)/(mean(x[(length(response_ids)+1):nrow(melanoma_PD1)])+0.01))
p_value=apply(ordered_GSVA,1,function(x) t.test(x[1:length(response_ids)],x[(length(response_ids)+1):nrow(melanoma_PD1)])$p.value)
result=as.data.frame(cbind(ordered_GSVA,response_Mean,non_response_Mean,FC,p_value))

cbind(rownames(result),result) %>%
  dplyr::filter(p_value<=0.01) -> sig_sets

rownames(sig_sets) = sig_sets[,1]
dplyr::select(sig_sets,2:(ncol(sig_sets)-4)) %>%
  apply(1,function(x) mean(x)) -> 


