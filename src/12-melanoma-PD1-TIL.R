library(magrittr)

#melanoma PD1 TIL
read.table("/data/miaoyr/immune_checkpoint/RNAseq_TIL_result/melanoma_pd1_fpkm_til_5_14",sep="\t",header = T) ->melanoma_PD1_TIL


readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Cancer=="melanoma")%>%
  dplyr::filter(Anti_target=="anti-PD1")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(Survival_time != "NA")%>%
  #dplyr::filter(SRA_Study == "SRP094781")%>%
  dplyr::select(Run,Response,Survival_time,Survival_status,Age,Gender) ->metadata_survival


tibble::rownames_to_column(melanoma_PD1_TIL)%>%
  merge(metadata_survival,.,by.x="Run",by.y="rowname")->Combined_data

# coxph(Surv(Survival_time, Survival_status)~InfiltrationScore, data = Combined_data)

for (j in 1:nrow(Combined_data)) {
  if(Combined_data$Survival_status[j]=="Dead"){
    Combined_data$Survival_status[j]="2"
  }else {
    Combined_data$Survival_status[j]="1"
  }
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

#Univariate Cox regression
covariates <- colnames(Combined_data)[-c(1,2,3,4)]

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
colnames(filter_res)[1]="TIL"
filter_res %>%
  dplyr::filter(p.value<=0.05)%>%
  dplyr::select(TIL)->TIL 
# 1   Tr1
# 2 nTreg

Multi_cox <- coxph(Surv(Survival_time, Survival_status) ~ Tr1 + nTreg, data =  Combined_data)
# > summary(Multi_cox)
# Call:
#   coxph(formula = Surv(Survival_time, Survival_status) ~ Tr1 + 
#           nTreg, data = Combined_data)
# 
# n= 77, number of events= 45 
# 
# coef exp(coef)  se(coef)      z Pr(>|z|)   
# Tr1    2.675045 14.513007  0.938159  2.851  0.00435 **
#   nTreg -6.312601  0.001813  2.193066 -2.878  0.00400 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Tr1   14.513007     0.0689 2.308e+00   91.2683
# nTreg  0.001813   551.4774 2.465e-05    0.1334
# 
# Concordance= 0.67  (se = 0.047 )
# Rsquare= 0.163   (max possible= 0.989 )
# Likelihood ratio test= 13.69  on 2 df,   p=0.001065
# Wald test            = 11.9  on 2 df,   p=0.002599
# Score (logrank) test = 11.95  on 2 df,   p=0.002543

Combined_data %>%
  dplyr::select(Run,Survival_time,Survival_status,Tr1,nTreg)->TIL_data
cutoff_Tr1=mean(TIL_data$Tr1)
cutoff_nTreg=mean(TIL_data$nTreg)

for (i in 1:nrow(TIL_data)) {
  if(TIL_data$Tr1[i]>=cutoff_Tr1){
    TIL_data$Tr1[i]="high"
  }else {
    TIL_data$Tr1[i]="low"
  }
  
  if(TIL_data$nTreg[i]>=cutoff_nTreg){
    TIL_data$nTreg[i]="high"
  }else {
    TIL_data$nTreg[i]="low"
  }
  
  
}


fit <- survfit(Surv(Survival_time, Survival_status) ~ Tr1, data = TIL_data)
pdf(file="/data/liull/immune-checkpoint-blockade/TIL/melanoma_PD1/Tr1.pdf")
ggsurvplot(fit, data = TIL_data, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")
dev.off()#0.0041

# fit <- survfit(Surv(Survival_time, Survival_status) ~ nTreg, data = TIL_data)
# pdf(file="/data/liull/immune-checkpoint-blockade/TIL/melanoma_PD1/nTreg.pdf")
# ggsurvplot(fit, data = TIL_data, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")
# dev.off()#0.057
write.table(TIL_data,"/data/liull/immune-checkpoint-blockade/TIL/melanoma_PD1/TIL_classify.txt",
            quote=FALSE,col.names = TRUE,row.names = FALSE,sep="\t")
