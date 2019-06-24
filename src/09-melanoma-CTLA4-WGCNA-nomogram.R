library(magrittr)
library(GSVA)
library(rms)
library(survival)

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(Run != "SRR3083584")%>%
  dplyr::select(Run,Second_Response_standard) ->metadata

metadata %>%
  dplyr::filter(Second_Response_standard %in% c("R","long-survival"))->R_metadata
metadata %>%
  dplyr::filter(Second_Response_standard %in% c("NR"))->NR_metadata

read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/melanoma_CTLA4_pretreatment_Symbol_FPKM_expr.txt",
           header = T,as.is = TRUE) ->melanoma_CTLA4_expr
melanoma_CTLA4_expr %>%
  dplyr::select(R_metadata$Run)->SRP011540_Symbol_FPKM_expr_R
melanoma_CTLA4_expr %>%
  dplyr::select(NR_metadata$Run)->SRP011540_Symbol_FPKM_expr_NR

write.table(SRP011540_Symbol_FPKM_expr_R,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4_second_standard/SRP011540_Symbol_FPKM_expr_R.txt",
            col.names = TRUE,row.names = TRUE,quote=FALSE,sep = "\t")
write.table(SRP011540_Symbol_FPKM_expr_NR,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_NR/melanoma_CTLA4_second_standard/SRP011540_Symbol_FPKM_expr_NR.txt",
            col.names = TRUE,row.names = TRUE,quote=FALSE,sep = "\t")

#WGCNA3.r SRP011540_Symbol_FPKM_expr_R.txt 0.90 10 0.25
#WGCNA3.r expr_NR.txt SRP011540_Symbol_FPKM_expr_NR.txt 10 0.25

read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4_second_standard/raw_module.assign.txt",header = T,as.is = TRUE) %>%
  dplyr::filter(module != "grey")->R_color_gene

R_list_name=unique(R_color_gene$module)
R_list_sets=list()

for (i in 1:length(R_list_name)) {
  
  dplyr::filter(R_color_gene,module == R_list_name[i])%>%
    dplyr::select(gene)%>%
    as.matrix()%>%
    as.character()%>%
    list()->R_list_sets[i]
  
}
names(R_list_sets)=paste("R_",R_list_name,sep="")

read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_NR/melanoma_CTLA4_second_standard/raw_module.assign.txt",header = T,as.is = TRUE) %>%
  dplyr::filter(module != "grey")->NR_color_gene

NR_list_name=unique(NR_color_gene$module)
NR_list_sets=list()

for (i in 1:length(NR_list_name)) {
  
  dplyr::filter(NR_color_gene,module == NR_list_name[i])%>%
    dplyr::select(gene)%>%
    as.matrix()%>%
    as.character()%>%
    list()->NR_list_sets[i]
  
}
names(NR_list_sets)=paste("NR_",NR_list_name,sep="")

list_sets=c(R_list_sets,NR_list_sets)

read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/melanoma_CTLA4_pretreatment_Symbol_FPKM_expr.txt",
           header = T,as.is = TRUE)->melanoma_CTLA4_expr

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::filter(Survival_time != "NA")%>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(Run != "SRR3083584")%>%
  dplyr::select(Run,Response,Survival_time,Survival_status,Age,Gender) ->metadata


ssgava_score <- gsva(as.matrix(melanoma_CTLA4_expr), list_sets, min.sz=1, max.sz=999999, method="ssgsea",kcdf="Gaussian")

t(ssgava_score) %>%
  as.data.frame(stringsAsFactors=FALSE)->t_ssgava_score
cbind(rownames(t_ssgava_score),t_ssgava_score)->t_ssgava_score
rownames(t_ssgava_score)=NULL
colnames(t_ssgava_score)[1]="Run"

Combined_data=merge(metadata[,c(1,2,3,4,5,6)],t_ssgava_score)

for (j in 1:nrow(Combined_data)) {
  if(Combined_data$Survival_status[j]=="dead"){
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
colnames(filter_res)[1]="module"
filter_res %>%
  dplyr::filter(p.value<=0.05)%>%
  dplyr::select(module)->selected_modules

# > selected_modules
# module
# 1     R_indianred3
# 2 R_paleturquoise4
# 3     R_powderblue
# 4     R_royalblue2
# 5     R_firebrick2
# 6  R_antiquewhite2
# 7        NR_purple
# 8 NR_mediumpurple3

Multi_cox <- coxph(Surv(Survival_time, Survival_status) ~ 
                     R_indianred3 + R_paleturquoise4 + R_powderblue + R_royalblue2 + R_firebrick2 +
                     R_antiquewhite2 + NR_purple + NR_mediumpurple3, data =  Combined_data)
# coef exp(coef)  se(coef)     z      p
# R_indianred3      2.37e+00  1.07e+01  2.77e+00  0.86 0.3916
# R_paleturquoise4  4.21e+00  6.77e+01  1.60e+00  2.63 0.0085
# R_powderblue      3.70e+00  4.03e+01  1.57e+00  2.36 0.0182
# R_royalblue2      4.18e+00  6.55e+01  1.63e+00  2.56 0.0103
# R_firebrick2     -4.44e+00  1.18e-02  2.85e+00 -1.56 0.1193
# R_antiquewhite2  -5.89e+00  2.78e-03  3.72e+00 -1.58 0.1140
# NR_purple        -2.25e+00  1.06e-01  6.84e+00 -0.33 0.7426
# NR_mediumpurple3 -1.14e+01  1.14e-05  4.12e+00 -2.76 0.0057
# 
# Likelihood ratio test=38.3  on 8 df, p=6.57e-06
# n= 41, number of events= 28 

#R_paleturquoise4 R_powderblue R_royalblue2 NR_mediumpurple3 plot---------------------------------------------------------------------
Combined_data %>%
  dplyr::select(Run,Survival_time,Survival_status,R_paleturquoise4)->Combined_data_plot

cutoff=mean(Combined_data_plot$R_paleturquoise4)
Combined_data_plot %>%
  dplyr::mutate(Class=if_else(R_paleturquoise4 <= cutoff , "lower","higher"))->Combined_data_plot
fit <- survfit(Surv(Survival_time, Survival_status) ~ Class, data = Combined_data_plot)

pdf(file="/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4_second_standard/R_paleturquoise4_survival.pdf")
ggsurvplot(fit, data = Combined_data_plot, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")
dev.off()
#--------------------------------------------------------------------------------------------------

write.table(list_sets$R_paleturquoise4,"/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4/nomogram/R_paleturquoise4_genes.txt",
            row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
write.table(list_sets$R_powderblue,"/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4/nomogram/R_powderblue_genes.txt",
            row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
write.table(list_sets$R_royalblue2,"/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4/nomogram/R_royalblue2_genes.txt",
            row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
write.table(list_sets$NR_mediumpurple3,"/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4/nomogram/NR_mediumpurple3_genes.txt",
            row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")

Combined_data %>%
  dplyr::select(Run,Response,Survival_time,Survival_status,R_paleturquoise4,R_powderblue,R_royalblue2,NR_mediumpurple3)->cox_selected
# write.table(SRP094781_cox_selected,
#             "/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram//SRP094781_cox_selected.txt",
#             col.names = TRUE,row.names = FALSE,quote=FALSE,sep = "\t")

cox_selected ->all_class
all_class%>%
  dplyr::filter(Response != "X")->all_class
ddist <- datadist(all_class)
options(datadist='ddist')

f_cph <- cph(Surv(Survival_time,Survival_status) ~
               Response +
               R_paleturquoise4 +
               R_powderblue +
               R_royalblue2 +
               NR_mediumpurple3, 
             data=all_class,surv=TRUE)
1-rcorrcens(Surv(Survival_time,Survival_status) ~ predict(f_cph), data =  all_class)#[[1]]
# Somers' Rank Correlation for Censored Data    Response variable:Surv(Survival_time, Survival_status)
# 
#                   C   Dxy  aDxy    SD     Z P   n
# predict(f_cph) 0.83 1.659 0.341 0.939 -9.76 1 -39

f_cph <- cph(Surv(Survival_time,Survival_status) ~
               Response +
               R_paleturquoise4 +
               R_powderblue +
               R_royalblue2 +
               NR_mediumpurple3, 
             data=all_class,surv=TRUE,x=TRUE, y=TRUE,time.inc=730)
surv <- Survival(f_cph)
nom_cph <- nomogram(f_cph, fun=list(function(x) surv(365, x),function(x) surv(730, x)),
                    fun.at = c(0.1,0.3,0.5,0.7,0.9),
                    funlabel=c("1-Year-Survival","2-Year-Survival"),
                    lp=F
)

jpeg('/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4_second_standard/nomogram/nomogram.jpeg',
     width = 1200, height = 800, units = "px", pointsize = 12,quality = 100,bg = "#e5ecff",res=100)
plot(nom_cph, xfrac=.3,cex.axis=.7,cex.var=0.9,lmgp=.05,tcl=0.25)
dev.off()

cal2 <- calibrate(f_cph, cmethod='KM', method="boot", u=730, m=12)

jpeg('/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4_second_standard/nomogram/calibration_curves.jpeg',
     width = 800, height = 800, units = "px", pointsize = 12,quality = 100,bg = "#e5ecff",res=100)
par(mar=c(8,5,3,2),cex = 1.0)
plot(cal2,lwd=2,lty=1,
     #errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0,1.0),ylim=c(0,1.0),
     xlab="Nomogram-Predicted Probability of 2-Year OS",
     ylab="Actual 2-Year OS",
     col=c(rgb(192,98,83,maxColorValue=255)),subtitles=FALSE)
abline(0,1,lty =3)
dev.off()
###

read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4_second_standard/nomogram/R_paleturquoise4_genes.txt",
           header = F)->R_paleturquoise4
read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4_second_standard/nomogram/R_powderblue_genes.txt",
           header = F)->R_powderblue
read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4_second_standard/nomogram/R_royalblue2_genes.txt",
           header = F)->R_royalblue2
read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/melanoma_CTLA4_second_standard/nomogram/NR_mediumpurple3_genes.txt",
           header = F)->NR_mediumpurple3
read.table("/data/liull/immune-checkpoint-blockade/TIL/melanoma_CTLA4/neoantigen_related_gene.txt",
           header = T)-> neoantigen_related_gene

#if devided by RECIST stabdard
# coef exp(coef)  se(coef)     z     p
# R_lightpink3      7.96e+00  2.85e+03  4.84e+00  1.65 0.100
# R_lavenderblush2  4.69e+00  1.08e+02  4.94e+00  0.95 0.343
# R_aliceblue       5.91e+00  3.68e+02  3.14e+00  1.88 0.060
# R_deeppink2       5.32e+00  2.05e+02  3.18e+00  1.68 0.094
# R_magenta4        1.53e+01  4.40e+06  9.81e+00  1.56 0.119
# NR_blue          -1.61e+01  1.02e-07  1.85e+01 -0.87 0.383
# NR_blue2          3.81e+00  4.53e+01  4.72e+00  0.81 0.419
# NR_darkorange     1.54e+00  4.65e+00  6.59e+00  0.23 0.815
# NR_indianred3    -9.27e+00  9.39e-05  9.95e+00 -0.93 0.352
# NR_plum2          4.89e+00  1.33e+02  3.20e+00  1.53 0.127
# NR_darkviolet    -1.09e+01  1.91e-05  6.65e+00 -1.63 0.102
# NR_maroon         4.16e+00  6.39e+01  5.59e+00  0.74 0.457
# 
# Likelihood ratio test=41.9  on 12 df, p=3.43e-05
# n= 41, number of events= 28 
