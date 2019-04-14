library(survminer)
library(survival)
#
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Survival_time != "NA")%>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%###
  dplyr::select(Run,Survival_time,Survival_status,Age,Gender) ->metadata#SRA all anti-PD1

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/PD1_removed_batch_expression.txt",header = T,as.is = TRUE)->expr  #removed batch effect expr

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
tibble::rownames_to_column(expr)%>%
  merge(relationship,.,by.x="Ensembl_ID",by.y="rowname")%>%
  dplyr::select(Symbol,metadata$Run)->expr2

factors=factor(expr2$Symbol)
merged_expression=tapply(expr2[,2],factors,median)
for (i in 3:ncol(expr2)) {
  temp=tapply(expr2[,i],factors,median)
  merged_expression=cbind(merged_expression,temp)
}
colnames(merged_expression)=colnames(expr2)[2:ncol(expr2)]


Collagen_deg=c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL9A3","COL12A1","COL15A1","COL16A1","CTSK","MMP2","MMP11")
ECM_deg=c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL9A3","COL12A1","COL15A1","COL16A1","CAPN12","CTSK","DCN","ELN","LAMC1","MMP2","MMP11","MMP17","NID1","HTRA1","BMP1","CAPN2","ADAMTS4","ADAMTS1")
ECM_org=c("CRTAP","FBLN5","P3H3","COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL9A3","COL12A1","COL15A1","COL16A1","VCAN","CAPN12","CTSK","DCN","ELN","FBLN1","TNC","ITGA6","ITGA5","ITGA7","AGRN","KDR","LAMA4","LAMC1","LOXL1","LTBP3","LUM","MFAP1","MFAP2","MFAP4","MMP2","MMP11","MMP14","MMP17","NID1","P4HB","PDGFB","PECAM1","PLEC","P3H2","HTRA1","PTPRS","P3H1","BMP1","SPARC","TGFB3","VWF","PXDN","MFAP5","COL18A1","CAPN2","COL27A1","SERPINH1","ADAM15","P4HA2","PLOD3","ADAMTS4","ADAMTS2","ADAMTS1")
Collagen_form=c("CRTAP","P3H3","COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL9A3","COL12A1","COL15A1","COL16A1","ITGA6","P4HB","PLEC","P3H1","BMP1","PXDN","COL27A1","SERPINH1","P4HA2","PLOD3","ADAMTS2")

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_ENSG.txt",header = T,as.is = TRUE)->gastric_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_ENSG.txt",header = T,as.is = TRUE)->gastric_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_up_ENSG.txt",header = T,as.is = TRUE)->melanoma_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_down_ENSG.txt",header = T,as.is = TRUE)->melanoma_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_up_ENSG.txt",header = T,as.is = TRUE)->melanoma_CTLA4_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_down_ENSG.txt",header = T,as.is = TRUE)->melanoma_CTLA4_down

intersect(gastric_PD1_down$Symbol,melanoma_PD1_down$Symbol)->PD1_down_melanoma_gastric


# down_melanoma_PD1_CTLA4_Interaction=c("UBE2E3","MRC2")
# down_melanoma_PD1_gastric=c("MMP2","PIR","CYBRD1","COL6A1","COL3A1")
# down_melanoma_CTLA4_gastric_PD1=c("NID1","ATP2B1","FSTL1","STXBP5","NDRG1")

union(Collagen_deg,ECM_deg) %>% union(ECM_org) %>% union(Collagen_form) %>% 
  union(PD1_down_melanoma_gastric) -> all_genes





#filter gene and sample ...
tibble::rownames_to_column(as.data.frame(merged_expression)) %>%
  dplyr::filter(rowname %in% all_genes)->expr3
rownames(expr3)=expr3$rowname
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

# > gene_symbol
# gene_symbol
# 1      CORO2B
# 2       MFAP2
# 3        NSG1
# 4        P3H1
# 5       P4HA2
# 6       RCOR2
# 7       RSPO4
# 8    SERPINH1

write.table(gene_symbol,"/data/liull/immune-checkpoint-blockade/survival/PD1_survival_symbol.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#Multivariate Cox regression analysis for significant gene_symbol
Multi_cox <- coxph(Surv(Survival_time, Survival_status) ~ COL16A1 + PIR + PDGFB + P3H3 +MFAP2 +TGFB3 + SERPINH1 + HTRA1 + UBE2E3 + CAPN12 + MMP17, data =  Combined_data)
#summary(res.cox)


#significant symbols' expression
Combined_data %>%
  dplyr::select(Run,as.character(gene_symbol$gene_symbol)) ->sig_gene_expr
sig_gene_expr=sig_gene_expr[,-1]
nrow(sig_gene_expr)
frame_to_plots=data.frame()
temp=NULL
for (i in 1:length(colnames(sig_gene_expr))) {
  cbind(sig_gene_expr[,i],rep(colnames(sig_gene_expr)[i],nrow(sig_gene_expr)))->temp
  rbind(frame_to_plots,temp)->frame_to_plots
}

frame_to_plots$V1=as.numeric(as.character(frame_to_plots$V1))
frame_to_plots$V2=as.character(frame_to_plots$V2)
colnames(frame_to_plots)=c("expression","Symbol")

ggboxplot(frame_to_plots, x="Symbol", y="expression", color = "Symbol")->PD1_symbl_expr
ggsave(
  filename = 'boxplot_PD1_symbl_expr.pdf',
  plot = PD1_symbl_expr,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/survival',
  width = 6,
  height = 6.8
)


#plot the survival plot
Combined_data %>%
  dplyr::select(Run,Survival_time,Survival_status,as.character(gene_symbol[,1]))->Combined_data_2

Mean=apply(Combined_data_2[,4:ncol(Combined_data_2)], 1,function(x) mean(x))
Combined_data_2 %>%
  dplyr::mutate(Mean=Mean) %>%
  dplyr::mutate(Class=rep("class",nrow(Combined_data_2)))%>%
  dplyr::select(Run,Survival_time,Survival_status,Mean,Class)->Combined_data_3

cutoff=mean(Mean)


for (i in 1:nrow(Combined_data_3)) {
  if(Combined_data_3$Mean[i]>=cutoff){
    Combined_data_3$Class[i]="high"
  }else {
    Combined_data_3$Class[i]="low"
  }
}

fit <- survfit(Surv(Survival_time, Survival_status) ~ Class, data = Combined_data_3)
pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/8_gene_survival.pdf")
ggsurvplot(fit, data = Combined_data_3, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")
dev.off()

#SERPINH1 survival plot
Combined_data %>%
  dplyr::select(Run,Survival_time,Survival_status,SERPINH1)->SERPINH1_data
cutoff=mean(SERPINH1_data$SERPINH1)

Combined_data_2 %>%
  dplyr::mutate(Mean=mean()) %>%
  dplyr::mutate(Class=rep("class",nrow(Combined_data_2)))%>%
  dplyr::select(Run,Survival_time,Survival_status,Mean,Class)->Combined_data_3




for (i in 1:nrow(SERPINH1_data)) {
  if(SERPINH1_data$SERPINH1[i]>=cutoff){
    SERPINH1_data$SERPINH1[i]="high"
  }else {
    SERPINH1_data$SERPINH1[i]="low"
  }
}


fit <- survfit(Surv(Survival_time, Survival_status) ~ SERPINH1, data = SERPINH1_data)
ggsurvplot(fit, data = SERPINH1_data, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")
