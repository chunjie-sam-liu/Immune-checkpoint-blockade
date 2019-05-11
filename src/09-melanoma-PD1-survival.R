library(survminer)
library(survival)


#Single Cox regression analysis for gene_symbol(ECM,collagen,intersect(melanoma_PD1_down,melanoma_CTLA4_down))----------------------------------
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Survival_time != "NA")%>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%###
  dplyr::select(Run,Response,Survival_time,Survival_status,Age,Gender) ->metadata#SRA all anti-PD1

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/PD1_removed_batch_expression.txt",header = T,as.is = TRUE)->expr  #removed batch effect expr

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
tibble::rownames_to_column(expr)%>%
  merge(relationship,.,by.x="Ensembl_ID",by.y="rowname")%>%
  dplyr::select(Symbol,metadata$Run)->expr2

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_log2CPM_expr.txt",
           sep="\t",header = T,as.is = TRUE)->merged_expression


Collagen_deg=c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL9A3","COL12A1","COL15A1","COL16A1","CTSK","MMP2","MMP11")
ECM_deg=c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL9A3","COL12A1","COL15A1","COL16A1","CAPN12","CTSK","DCN","ELN","LAMC1","MMP2","MMP11","MMP17","NID1","HTRA1","BMP1","CAPN2","ADAMTS4","ADAMTS1")
ECM_org=c("CRTAP","FBLN5","P3H3","COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL9A3","COL12A1","COL15A1","COL16A1","VCAN","CAPN12","CTSK","DCN","ELN","FBLN1","TNC","ITGA6","ITGA5","ITGA7","AGRN","KDR","LAMA4","LAMC1","LOXL1","LTBP3","LUM","MFAP1","MFAP2","MFAP4","MMP2","MMP11","MMP14","MMP17","NID1","P4HB","PDGFB","PECAM1","PLEC","P3H2","HTRA1","PTPRS","P3H1","BMP1","SPARC","TGFB3","VWF","PXDN","MFAP5","COL18A1","CAPN2","COL27A1","SERPINH1","ADAM15","P4HA2","PLOD3","ADAMTS4","ADAMTS2","ADAMTS1")
Collagen_form=c("CRTAP","P3H3","COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL9A3","COL12A1","COL15A1","COL16A1","ITGA6","P4HB","PLEC","P3H1","BMP1","PXDN","COL27A1","SERPINH1","P4HA2","PLOD3","ADAMTS2")

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_ENSG.txt",header = T,as.is = TRUE)->gastric_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_ENSG.txt",header = T,as.is = TRUE)->gastric_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_up_ENSG.txt",header = T,as.is = TRUE)->melanoma_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_down_ENSG.txt",header = T,as.is = TRUE)->melanoma_PD1_down
# read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_up_ENSG.txt",header = T,as.is = TRUE)->melanoma_CTLA4_up
# read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_down_ENSG.txt",header = T,as.is = TRUE)->melanoma_CTLA4_down

intersect(gastric_PD1_down$Symbol,melanoma_PD1_down$Symbol)->PD1_down_melanoma_gastric

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

write.table(gene_symbol,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/PD1_survival_symbol.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)

#Multivariate Cox regression analysis for significant gene_symbol-----------------------------------------------------------------------------------
Multi_cox <- coxph(Surv(Survival_time, Survival_status) ~ CORO2B + MFAP2 + NSG1 + P3H1 +P4HA2 +RCOR2 + RSPO4 + SERPINH1, data =  Combined_data)
# > summary(Multi_cox)
# Call:
#   coxph(formula = Surv(Survival_time, Survival_status) ~ CORO2B + 
#           MFAP2 + NSG1 + P3H1 + P4HA2 + RCOR2 + RSPO4 + SERPINH1, data = Combined_data)
# 
# n= 77, number of events= 45 
# 
# coef exp(coef)  se(coef)     z Pr(>|z|)  
# CORO2B   0.0897150 1.0938625 0.0795505 1.128   0.2594  
# MFAP2    0.1529204 1.1652323 0.1082760 1.412   0.1579  
# NSG1     0.1499755 1.1618058 0.0828150 1.811   0.0701 .
# P3H1     0.0007536 1.0007539 0.2139418 0.004   0.9972  
# P4HA2    0.1367937 1.1465915 0.1460346 0.937   0.3489  
# RCOR2    0.0494211 1.0506627 0.1082023 0.457   0.6479  
# RSPO4    0.0565511 1.0581807 0.0554413 1.020   0.3077  
# SERPINH1 0.0495130 1.0507592 0.1763723 0.281   0.7789  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# CORO2B       1.094     0.9142    0.9359     1.278
# MFAP2        1.165     0.8582    0.9424     1.441
# NSG1         1.162     0.8607    0.9877     1.367
# P3H1         1.001     0.9992    0.6580     1.522
# P4HA2        1.147     0.8722    0.8612     1.527
# RCOR2        1.051     0.9518    0.8499     1.299
# RSPO4        1.058     0.9450    0.9492     1.180
# SERPINH1     1.051     0.9517    0.7437     1.485
# 
# Concordance= 0.672  (se = 0.047 )
# Rsquare= 0.238   (max possible= 0.989 )
# Likelihood ratio test= 20.91  on 8 df,   p=0.00738
# Wald test            = 16.76  on 8 df,   p=0.0327
# Score (logrank) test = 17.93  on 8 df,   p=0.02176


#significant symbols' expression---------------------------------------------------------------------------------------
Combined_data %>%
  dplyr::select(Run,as.character(gene_symbol$gene_symbol)) ->sig_gene_expr
merge(metadata,sig_gene_expr)%>%
  dplyr::select(-Survival_time,-Survival_status,-Age,-Gender)%>%
  dplyr::filter(Response != "NE")->sig_gene_expr
sig_gene_expr$Response%>%
  gsub("^PD$", "NR",. )%>%
  gsub("^SD$", "NR", .)%>%
  gsub("^PR$", "R", .)%>%
  gsub("^CR$", "R", .)->sig_gene_expr$Response

sig_gene_expr=sig_gene_expr[,-1]
nrow(sig_gene_expr)
frame_to_plots=data.frame()
temp=NULL
for (i in 2:length(colnames(sig_gene_expr))) {
  cbind(sig_gene_expr[,i],rep(colnames(sig_gene_expr)[i],nrow(sig_gene_expr)))%>%
    cbind(sig_gene_expr[,1])->temp
  rbind(frame_to_plots,temp)->frame_to_plots
}

frame_to_plots$V1=as.numeric(as.character(frame_to_plots$V1))
frame_to_plots$V2=as.character(frame_to_plots$V2)
frame_to_plots$V3=as.character(frame_to_plots$V3)
colnames(frame_to_plots)=c("expression","Symbol","Response")

ggboxplot(frame_to_plots, x="Symbol", y="expression", color = "Response")->PD1_symbl_expr
ggsave(
  filename = 'boxplot_PD1_symbl_expr.pdf',
  plot = PD1_symbl_expr,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/',
  width = 10,
  height = 6.8
)


#plot the survival plot------------------------------------------------------------------------------------------------------
Combined_data %>%
  dplyr::select(Run,Survival_time,Survival_status,as.character(gene_symbol[,1]))->Combined_data_2

Mean=apply(Combined_data_2[,4:ncol(Combined_data_2)],1,function(x) mean(x))
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

Combined_data_3[,c(1,2,5)] %>%
  merge(metadata)->eight_gene_class
colnames(eight_gene_class)[3]="eight_gene_class"
write.table(eight_gene_class,
            "/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/8_gene_class.txt",
            col.names = TRUE,row.names = FALSE,quote=FALSE)

fit <- survfit(Surv(Survival_time, Survival_status) ~ Class, data = Combined_data_3)
pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/8_gene_survival.pdf")
ggsurvplot(fit, data = Combined_data_3, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")
dev.off()

#single gene (SERPINH1) survival plot-------------------------------------------------------------------------------------
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


#heatmap------------------------------------------------------------------------------------------------------------------
dplyr::filter(metadata,Response %in% c("CR","PR")) -> response#
dplyr::filter(metadata,Response %in% c("SD","PD")) -> non_response#


tibble::rownames_to_column(expr3) %>%
  dplyr::select(rowname,response$Run,non_response$Run)%>%
  dplyr::filter(rowname %in% gene_symbol$gene_symbol)->filtered_gene_expr
rownames(filtered_gene_expr)=filtered_gene_expr$rowname
filtered_gene_expr[,-1]->filtered_gene_expr

apply(filtered_gene_expr, 1, scale) ->scaled_expr
rownames(scaled_expr)=colnames(filtered_gene_expr)
scaled_expr=t(scaled_expr)


df = data.frame(type = c(rep("response", nrow(response)), rep("non_response", nrow(non_response))))
ha = HeatmapAnnotation(df = df,col = list(type = c("response" =  "tomato", "non_response" = "steelblue")))

pdf(file="/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/8_gene_heatmap.pdf")
Heatmap(scaled_expr,name="Color_key",top_annotation = ha,cluster_columns = FALSE,column_names_gp = gpar(fontsize = 2),row_names_gp = gpar(fontsize = 6),col=colorRamp2(c(-2, 0, 2), c("green", "black", "red")))
dev.off()
