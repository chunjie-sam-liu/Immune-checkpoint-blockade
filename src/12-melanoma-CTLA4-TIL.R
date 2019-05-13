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