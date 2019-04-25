library(magrittr)

#filter the TF associated with survival with Univariate Cox regression------------------------------------------------------------------------------------------
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

#TF-target in modules----------------------------------------------------------------------------------------------------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/TF_symbol.txt",header = T,as.is = TRUE)->TF_symbol
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/NR_white.txt",header = F,as.is = TRUE)->NR_white#0 interaction
colnames(NR_white)="module_genes"
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/NR_orangered4.txt",header = F,as.is = TRUE)->NR_orangered4#0 TF
colnames(NR_orangered4)="module_genes"
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/NR_midnightblue.txt",header = F,as.is = TRUE)->NR_midnightblue#*****
colnames(NR_midnightblue)="module_genes"
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/R_skyblue3.txt",header = F,as.is = TRUE)->R_skyblue3#0 interaction
colnames(R_skyblue3)="module_genes"


read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/chip_tf_to_coding_gene.txt",header = T,sep="\t",as.is = TRUE)->TF_targets

dplyr::filter(TF_symbol,gene_symbol %in% NR_midnightblue$module_genes)->TFs_in_module
dim(TFs_in_module)
dplyr::mutate(NR_midnightblue,target=rep("NA",nrow(NR_midnightblue)))->NR_midnightblue

Targets_relationship=data.frame()
colnames(Targets_relationship)=c("TF","Targets")
for (i in 1:nrow(TFs_in_module)) {
  
  every_TF=TFs_in_module$gene_symbol[i]
  print(every_TF)
  dplyr::filter(TF_targets,TF == every_TF)%>%
    dplyr::select(target)%>%
    as.matrix()%>%
    as.character()->a
  dplyr::filter(TF_targets,target == every_TF)%>%
    dplyr::select(TF)%>%
    as.matrix()%>%
    as.character()->b
  union(a,b)%>%unique()->every_TF_targets
  intersect(every_TF_targets,NR_midnightblue$module_genes)->targets_in_module
  data.frame(TF=rep(every_TF,length(targets_in_module)),Targets=targets_in_module) %>%
    rbind(Targets_relationship) ->Targets_relationship
  # intersect(every_TF_targets,NR_midnightblue$module_genes)->m
  # print(m)
  
}

# [1] 3
# [1] 9
# [1] 1
# [1] 2

# > Targets_relationship
# TF Targets
# 1  ZNF623   RAD21
# 2  ZNF623  ZBTB10
# 3    ZHX1   RAD21
# 4  ZBTB10  VCPIP1
# 5  ZBTB10   RMDN1
# 6  ZBTB10  TMEM67
# 7  ZBTB10  RAD54B
# 8  ZBTB10  VPS13B
# 9  ZBTB10   YWHAZ
# 10 ZBTB10  TMEM65
# 11 ZBTB10  ZNF623
# 12 ZBTB10   RAD21
# 13  TERF1   MTFR1
# 14  TERF1  VPS13B
# 15  TERF1   RAD21













