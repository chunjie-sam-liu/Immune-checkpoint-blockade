library(survminer)
library(survival)
#

read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_CTLA4_removed_batch_expression.txt",header = T,as.is = TRUE)->expr  # CTLA4
cbind(rownames(expr),expr)->expr
colnames(expr)[1]="ensembl_ID"
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,expr,by.x="Ensembl_ID",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% expr$ensembl_ID) %>%
  dplyr::select(-Ensembl_ID,-GeneID)->expr2

Collagen_deg=c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL9A3","COL12A1","COL15A1","COL16A1","CTSK","MMP2","MMP11")
ECM_deg=c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL9A3","COL12A1","COL15A1","COL16A1","CAPN12","CTSK","DCN","ELN","LAMC1","MMP2","MMP11","MMP17","NID1","HTRA1","BMP1","CAPN2","ADAMTS4","ADAMTS1")
ECM_org=c("CRTAP","FBLN5","P3H3","COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL9A3","COL12A1","COL15A1","COL16A1","VCAN","CAPN12","CTSK","DCN","ELN","FBLN1","TNC","ITGA6","ITGA5","ITGA7","AGRN","KDR","LAMA4","LAMC1","LOXL1","LTBP3","LUM","MFAP1","MFAP2","MFAP4","MMP2","MMP11","MMP14","MMP17","NID1","P4HB","PDGFB","PECAM1","PLEC","P3H2","HTRA1","PTPRS","P3H1","BMP1","SPARC","TGFB3","VWF","PXDN","MFAP5","COL18A1","CAPN2","COL27A1","SERPINH1","ADAM15","P4HA2","PLOD3","ADAMTS4","ADAMTS2","ADAMTS1")
Collagen_form=c("CRTAP","P3H3","COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL9A3","COL12A1","COL15A1","COL16A1","ITGA6","P4HB","PLEC","P3H1","BMP1","PXDN","COL27A1","SERPINH1","P4HA2","PLOD3","ADAMTS2")


down_melanoma_PD1_CTLA4_Interaction=c("UBE2E3","MRC2")
down_melanoma_PD1_gastric=c("MMP2","PIR","CYBRD1","COL6A1","COL3A1")
down_melanoma_CTLA4_gastric_PD1=c("NID1","ATP2B1","FSTL1","STXBP5","NDRG1")

union(Collagen_deg,ECM_deg) %>% union(ECM_org) %>% union(Collagen_form) %>% 
  union(down_melanoma_PD1_CTLA4_Interaction)%>%union(down_melanoma_PD1_gastric)%>%
  union(down_melanoma_CTLA4_gastric_PD1) -> all_genes



readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Survival_time != "NA")%>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%###
  dplyr::select(Run,Survival_time,Survival_status,Age,Gender) ->metadata#SRA all anti-PD1,dnGAP all anti-CTLA4

#filter gene and sample ...
expr2 %>%
  dplyr::filter(Symbol %in% all_genes)%>%
  dplyr::select(Symbol,metadata$Run)->expr3
rownames(expr3)=expr3$Symbol
expr3=expr3[,-1]

expr4=as.data.frame(t(expr3))
cbind(rownames(expr4),expr4)->all_sets
rownames(all_sets)=NULL
colnames(all_sets)[1]="Run"

Combined_data=merge(metadata,all_sets)

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
  dplyr::filter(p.value<=0.05)%>%
  dplyr::select(gene_symbol)->gene_symbol
write.table(gene_symbol,"/data/liull/immune-checkpoint-blockade/survival/CTLA4_survival_symbol.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#Multivariate Cox regression analysis for significant gene_symbol
Multi_cox <- coxph(Surv(Survival_time, Survival_status) ~ COL16A1 + PIR + PDGFB + P3H3 +MFAP2 +TGFB3 + SERPINH1 + HTRA1 + UBE2E3 + CAPN12 + MMP17, data =  Combined_data)