read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/Homo_sapiens_TF",header = T,as.is = TRUE,sep="\t")%>%
  dplyr::select(Symbol)->human_TF

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/pre_PD1_filtered_symbol_expr.txt",header = T,as.is = TRUE)->pre_PD1_expr
tibble::rownames_to_column(pre_PD1_expr)%>%
  dplyr::filter(rowname %in% human_TF$Symbol)->TF_expr
rownames(TF_expr)=TF_expr$rowname
TF_expr=TF_expr[,-1]
t_TF_expr=as.data.frame(t(TF_expr))
cbind(rownames(t_TF_expr),t_TF_expr)->t_TF_expr
rownames(t_TF_expr)=NULL
colnames(t_TF_expr)[1]="Run"

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Survival_time != "NA")%>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%###
  dplyr::select(Run,Response,Survival_time,Survival_status,Age,Gender) ->metadata

Combined_data=merge(metadata[,c(1,3,4,5,6)],t_TF_expr)

for (j in 1:nrow(Combined_data)) {
  if(Combined_data$Survival_status[j]=="Dead"){
    Combined_data$Survival_status[j]="2"
  }else {
    Combined_data$Survival_status[j]="1"
  }
}

for (j in 1:nrow(Combined_data)) {
  if(Combined_data$Gender[j]=="female"){
    Combined_data$Gender[j]="2"
  }else if(Combined_data$Gender[j]=="male") {
    Combined_data$Gender[j]="1"
  }
}

Combined_data$Survival_time=as.numeric(Combined_data$Survival_time)
Combined_data$Survival_status=as.numeric(Combined_data$Survival_status)
Combined_data$Age=as.numeric(Combined_data$Age)
Combined_data$Gender=as.numeric(Combined_data$Gender)



gsub("-", "_",colnames(Combined_data))->colnames(Combined_data)
  
#Univariate Cox regression
covariates <- colnames(Combined_data)[-c(1,2,3)]

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(Survival_time, Survival_status)~',x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = Combined_data)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })


res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)),stringsAsFactors=FALSE)
res$p.value=as.numeric(res$p.value)

cbind(rownames(res),res)->filter_res
colnames(filter_res)[1]="gene_symbol"
filter_res %>%
  dplyr::filter(p.value<=0.2)%>%
  dplyr::select(gene_symbol)->TF_symbol#285
write.table(TF_symbol,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/TF_symbol.txt",quote = FALSE,col.names = TRUE,row.names = FALSE)
