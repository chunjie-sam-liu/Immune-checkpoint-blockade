#TIDE predict response

read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/SRP070710_Symbol_FPKM_expr.txt",
           header = T,as.is = TRUE) ->SRP070710_FPKM_expr

for (i in 1:ncol(SRP070710_FPKM_expr)) {
  
  SRP070710_FPKM_expr[,i][is.na(SRP070710_FPKM_expr[,i])]=0
  log2(SRP070710_FPKM_expr[,i]+1)->SRP070710_FPKM_expr[,i]

}

normalize.quantiles(as.matrix(SRP070710_FPKM_expr))->quantiled_expr
rownames(quantiled_expr)=rownames(SRP070710_FPKM_expr)
colnames(quantiled_expr)=colnames(SRP070710_FPKM_expr)
quantiled_expr=as.data.frame(quantiled_expr,stringsAsFactors = FALSE)

Means=apply(quantiled_expr, 1, function(x) mean(x))
Subtracted_expr = data.frame(row.names = rownames(quantiled_expr))
for (i in 1:ncol(quantiled_expr)) {
  
  quantiled_expr[,i]-Means -> Subtracted_expr[,i]
  
}
colnames(Subtracted_expr)=colnames(quantiled_expr)

write.table(Subtracted_expr,"/data/liull/immune-checkpoint-blockade/TIDE/Subtracted_SRP070710.txt",
            row.names = TRUE,col.names = TRUE,quote = FALSE,sep="\t")




readr::read_csv("/data/liull/immune-checkpoint-blockade/TIDE/SRP070710_tide_result.csv")->tide_result
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(SRA_Study=="SRP070710") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(Run,Response) ->metadata

metadata$Response %>%
  gsub("CR","R",.)%>%
  gsub("PR","R",.)%>%
  gsub("SD","NR",.)%>%
  gsub("PD","NR",.)->metadata$Response

tide_result %>%
  dplyr::select(Patient,TIDE)%>%
  merge(metadata,by.x = "Patient",by.y = "Run")->AUC_matrix


roc(AUC_matrix$Response, AUC_matrix$TIDE, plot=TRUE, print.auc=TRUE)#0.7444
