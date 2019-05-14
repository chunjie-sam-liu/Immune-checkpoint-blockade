library(magrittr)

#melanoma PD1 TIL
read.table("/data/miaoyr/immune_checkpoint/RNAseq_TIL_result/melanoma_CTLA4_fpkm_til",sep="\t",header = T) ->melanoma_CTLA4_TIL

#survival
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP") %>%
  dplyr::filter(Cancer=="melanoma")%>%
  dplyr::filter(Anti_target=="anti-CTLA4")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(Survival_time != "NA")%>%
  dplyr::filter(Run != "SRR3083584") %>%
  dplyr::select(Run,Response,Survival_time,Survival_status,Age,Gender) ->metadata_survival


tibble::rownames_to_column(melanoma_CTLA4_TIL) %>%
  merge(metadata_survival,.,by.x="Run",by.y="rowname")->Combined_data

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
# coxph(Surv(Survival_time, Survival_status)~Neutrophil, data = Combined_data) p=0.2

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
# 1         Cytotoxic
# 2              Th17
# 3        Macrophage
# 4                NK
# 5        Neutrophil
# 6             CD8_T
# 7 InfiltrationScore  #######count

# 1   Tr1
# 2 nTreg
# 3  MAIT   ########FPKM
Multi_cox <- coxph(Surv(Survival_time, Survival_status) ~ Tr1 + nTreg + MAIT, data =  Combined_data)
# > summary(Multi_cox)
# Call:
#   coxph(formula = Surv(Survival_time, Survival_status) ~ Tr1 + 
#           nTreg + MAIT, data = Combined_data)
# 
# n= 41, number of events= 41 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)  
# Tr1    2.32411  10.21759  1.22253  1.901   0.0573 .
# nTreg -3.78874   0.02262  2.49368 -1.519   0.1287  
# MAIT  -2.42008   0.08891  1.30663 -1.852   0.0640 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# Tr1    10.21759    0.09787 0.9305300   112.193
# nTreg   0.02262   44.20078 0.0001706     3.000
# MAIT    0.08891   11.24673 0.0068670     1.151
# 
# Concordance= 0.716  (se = 0.054 )
# Rsquare= 0.27   (max possible= 0.996 )
# Likelihood ratio test= 12.91  on 3 df,   p=0.004841
# Wald test            = 12.92  on 3 df,   p=0.004824
# Score (logrank) test = 13.46  on 3 df,   p=0.003747

Combined_data %>%
  dplyr::select(Run,Survival_time,Survival_status,Tr1,nTreg,MAIT)->TIL_data
cutoff_Tr1=mean(TIL_data$Tr1)
cutoff_nTreg=mean(TIL_data$nTreg)
cutoff_MAIT=mean(TIL_data$MAIT)

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

  if(TIL_data$MAIT[i]>=cutoff_MAIT){
    TIL_data$MAIT[i]="high"
  }else {
    TIL_data$MAIT[i]="low"
  }
  
  
}


fit <- survfit(Surv(Survival_time, Survival_status) ~ Tr1, data = TIL_data)
pdf(file="/data/liull/immune-checkpoint-blockade/TIL/melanoma_CTLA4/Tr1.pdf")
ggsurvplot(fit, data = TIL_data, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")
dev.off()#0.0009

# fit <- survfit(Surv(Survival_time, Survival_status) ~ nTreg, data = TIL_data)
# pdf(file="/data/liull/immune-checkpoint-blockade/TIL/melanoma_CTLA4/nTreg.pdf")
# ggsurvplot(fit, data = TIL_data, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")
# dev.off()# 0.095

fit <- survfit(Surv(Survival_time, Survival_status) ~ MAIT, data = TIL_data)
pdf(file="/data/liull/immune-checkpoint-blockade/TIL/melanoma_CTLA4/MAIT.pdf")
ggsurvplot(fit, data = TIL_data, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")
dev.off() #0.048

