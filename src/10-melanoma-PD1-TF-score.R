library(magrittr)

#filter the TF associated with survival with Univariate Cox regression------------------------------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/Homo_sapiens_TF",header = T,as.is = TRUE,sep="\t")%>%
  dplyr::select(Symbol)->human_TF

read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/SRP094781_Symbol_FPKM_expr.txt",header = T,as.is = TRUE)->SRP094781_expr
tibble::rownames_to_column(SRP094781_expr)%>%
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
  dplyr::filter(SRA_Study == "SRP094781") %>%
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
covariates <- colnames(Combined_data)[-c(1,2,3,4,5)]

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


dim(TFs_in_black)
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








#new
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/Homo_sapiens_TF",
           header = T,as.is = TRUE,sep="\t")->all_TF
read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/NR_black_genes.txt",
           header = F,as.is = TRUE,sep="\t")->black_genes
read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/NR_ivory_genes.txt",
           header = F,as.is = TRUE,sep="\t")->ivory_genes

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/chip_tf_to_coding_gene.txt",header = T,sep="\t",as.is = TRUE)->TF_targets

dplyr::filter(all_TF,Symbol %in% black_genes$V1)->TFs_in_black

TF_targets %>%
  dplyr::filter(TF %in% TFs_in_black$Symbol) %>%
  dplyr::select(TF,target)->black_TF_targets_1

TF_targets %>%
  dplyr::filter(target %in% TFs_in_black$Symbol) %>%
  dplyr::select(target,TF)->black_TF_targets_2
colnames(black_TF_targets_2)=c("TF","target")
rbind(black_TF_targets_1,black_TF_targets_2) %>%
  unique()->black_TF_targets


# for (i in 1:nrow(TFs_in_black)) {
#   black_TF_targets %>%
#     dplyr::filter(TF == TFs_in_black$Symbol[i])->one_TF_targets
#   # if(nrow(one_TF_targets)==0){
#   #   black_TF_targets %>%
#   #     dplyr::filter(target == TFs_in_black$Symbol[i])->one_TF_targets
#   # }
#   
#   Num=length(intersect(one_TF_targets$target,black_genes$V1))
#   print(Num)
#   }
  
  


for (i in 1:nrow(TFs_in_black)) {
  black_TF_targets %>%
    dplyr::filter(TF == TFs_in_black$Symbol[i])->one_TF_targets
  if(nrow(one_TF_targets)==0){
    black_TF_targets %>%
      dplyr::filter(target == TFs_in_black$Symbol[i])->one_TF_targets
  }
  
  Num=length(intersect(one_TF_targets$target,black_genes$V1))
  if(Num > 90){
    TFs_in_black$Symbol[i] -> TF_signif
    write.table(intersect(one_TF_targets$target,black_genes$V1),
                paste("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/black_",TF_signif,".txt",sep = ""),
                row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
  }
  
  
}



dplyr::filter(all_TF,Symbol %in% ivory_genes$V1)->TFs_in_ivory

TF_targets %>%
  dplyr::filter(target %in% TFs_in_ivory$Symbol) %>%
  dplyr::select(TF,target)->ivory_TF_targets

for (i in 1:nrow(TFs_in_ivory)) {
  ivory_TF_targets %>%
    dplyr::filter(target == TFs_in_ivory$Symbol[i])->one_TF_targets
  if(nrow(one_TF_targets)==0){
    ivory_TF_targets %>%
      dplyr::filter(TF == TFs_in_ivory$Symbol[i])->one_TF_targets
  }
  
  Num=length(intersect(one_TF_targets$TF,black_genes$V1))
  if(Num > 90){
    TFs_in_ivory$Symbol[i] -> TF_signif
    write.table(intersect(one_TF_targets$target,black_genes$V1),
                paste("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/black_",TF_signif,".txt",sep = ""),
                row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
  }
  
  
}

# 1,5,1,4


read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/black_BHLHE40.txt",
           header = F,as.is = TRUE,sep="\t")->black_BHLHE40
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/black_RCOR1.txt",
           header = F,as.is = TRUE,sep="\t")->black_RCOR1
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/black_SMAD4.txt",
           header = F,as.is = TRUE,sep="\t")->black_SMAD4
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/black_TFAP2C.txt",
           header = F,as.is = TRUE,sep="\t")->black_TFAP2C
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/score_with_TF/black_TFDP1.txt",
           header = F,as.is = TRUE,sep="\t")->black_TFDP1

union(black_BHLHE40$V1,black_RCOR1$V1) %>%
  union(black_SMAD4$V1) %>% 
  union(black_TFAP2C$V1) %>%
  union(black_TFDP1$V1) %>%
  union(c("BHLHE40","RCOR1","SMAD4","TFAP2C","TFDP1"))->all_targets
list(all_targets) ->TF_targets_list[1]
names(TF_targets_list)="balck_TF"


