#WGCNA3.r expr_R.txt 0.85 30 0.25
#WGCNA3.r expr_NR.txt 0.85 30 0.25
#with the sets WGCNA get,do Univariate Cox regression----------------------------------------------------------------------------------
library(magrittr)
library(GSVA)
library(rms)


read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/raw_module.assign.txt",header = T,as.is = TRUE) %>%
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

read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_NR/SRP094781/raw_module.assign.txt",header = T,as.is = TRUE) %>%
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


read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/SRP094781_Symbol_FPKM_expr.txt",header = T,as.is = TRUE)->SRP094781_expr

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(SRA_Study == "SRP094781") %>%
  dplyr::filter(Survival_time != "NA")%>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(Run,Response,Survival_time,Survival_status,Age,Gender) ->metadata


ssgava_score <- gsva(as.matrix(SRP094781_expr), list_sets, min.sz=1, max.sz=999999, method="ssgsea",kcdf="Gaussian")

t(ssgava_score) %>%
  as.data.frame(stringsAsFactors=FALSE)->t_ssgava_score
cbind(rownames(t_ssgava_score),t_ssgava_score)->t_ssgava_score
rownames(t_ssgava_score)=NULL
colnames(t_ssgava_score)[1]="Run"

Combined_data=merge(metadata[,c(1,2,3,4)],t_ssgava_score)

for (j in 1:nrow(Combined_data)) {
  if(Combined_data$Survival_status[j]=="Dead"){
    Combined_data$Survival_status[j]="2"
  }else {
    Combined_data$Survival_status[j]="1"
  }
}

Combined_data$Survival_time=as.numeric(Combined_data$Survival_time)
Combined_data$Survival_status=as.numeric(Combined_data$Survival_status)

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
# module
# 1          R_yellow *
# 2       R_turquoise
# 3      R_darkorange 
# 4           R_white *
# 5          NR_black *
# 6 NR_palevioletred3 *
# 7          NR_ivory *
Multi_cox <- coxph(Surv(Survival_time, Survival_status) ~ R_yellow + R_turquoise + R_darkorange + R_white + NR_black + NR_palevioletred3 + NR_ivory, data =  Combined_data)
# coef exp(coef)  se(coef)     z       p
# R_yellow           1.05e+00  2.85e+00  6.18e+00  0.17 0.86562
# R_turquoise        6.71e+00  8.18e+02  3.93e+00  1.71 0.08753
# R_darkorange      -1.23e+00  2.93e-01  2.63e+00 -0.47 0.64083
# R_white            1.11e+01  6.61e+04  5.57e+00  1.99 0.04643
# NR_black          -1.72e+01  3.50e-08  6.07e+00 -2.83 0.00466
# NR_palevioletred3 -8.34e-02  9.20e-01  2.57e+00 -0.03 0.97411
# NR_ivory          -1.37e+01  1.16e-06  3.81e+00 -3.59 0.00033


#NR_black NR_ivory R_white

Combined_data %>%
  dplyr::select(Run,Response,Survival_time,Survival_status,NR_black,NR_ivory,R_white)->SRP094781_cox_selected
write.table(SRP094781_cox_selected,
            "/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram//SRP094781_cox_selected.txt",
            col.names = TRUE,row.names = FALSE,quote=FALSE,sep = "\t")

SRP094781_cox_selected ->all_class
all_class%>%
  dplyr::filter(Response != "NE")->all_class
ddist <- datadist(all_class)
options(datadist='ddist')
modules=colnames(all_class)[-c(1,2,3,4)]



#Survival_time,Survival_status ~ modules---------------------------------------------------
# for (i in 1:nrow(all_class)) {
#   if(all_class$Survival_status[i]=="Dead"){
#     all_class$Survival_status[i]=2
#   }else {
#     all_class$Survival_status[i]=1
#   }
# }
all_class$Survival_status=as.numeric(all_class$Survival_status)
all_class$Survival_time=as.numeric(all_class$Survival_time)

# #psm(up order)
# f_psm <- psm(Surv(Survival_time,Survival_status) ~ 
#                Response +
#                NR_black +
#                NR_ivory, 
#              data=all_class, x=T, y=T, dist='lognormal',time.inc=365)
# rcorrcens(Surv(Survival_time,Survival_status) ~ predict(f_psm), data =  all_class)
# # Somers' Rank Correlation for Censored Data    Response variable:Surv(Survival_time, Survival_status)
# # 
# # C   Dxy  aDxy    SD    Z P  n
# # predict(f_psm) 0.809 0.618 0.618 0.065 9.46 0 49
# 
# med  <- Quantile(f_psm)
# surv <- Survival(f_psm)  # This would also work if f was from cph
# nom_psm <- nomogram(f_psm, fun=list(function(x) surv(365, x),
#                                     function(x) surv(1095, x)),
#                     fun.at = c(0.1,0.2,0.4,0.6,0.8,0.9,1),
#                     funlabel=c("1-Year Survival Probability", 
#                                "3-Year Survival Probability"))
# plot(nom_psm, xfrac=.4,cex.axis=.7)
# 
# 
# cal1 <- calibrate(f_psm, cmethod='KM', method="boot", u=365, m=10, B=228)
# par(mar=c(8,5,3,2),cex = 1.0)
# plot(cal1,lwd=2,lty=1,
#      #errbar.col=c(rgb(0,118,192,maxColorValue=255)),
#      #xlim=c(0.25,0.6),ylim=c(0.15,0.70),
#      xlab="Nomogram-Predicted Probability of 1-Year DFS",
#      ylab="Actual 1-Year DFS (proportion)",
#      col=c(rgb(192,98,83,maxColorValue=255)))


#cph(down order)
f_cph <- cph(Surv(Survival_time,Survival_status) ~
               Response +
               NR_black +
               NR_ivory, 
             data=all_class,surv=TRUE)
1-rcorrcens(Surv(Survival_time,Survival_status) ~ predict(f_cph), data =  all_class)#[[1]]

# Somers' Rank Correlation for Censored Data    Response variable:Surv(Survival_time, Survival_status)
# 
#                   C  Dxy aDxy    SD     Z P   n
# predict(f_cph) 0.81 1.62 0.38 0.932 -8.17 1 -48

f_cph <- cph(Surv(Survival_time,Survival_status) ~
               Response +
               NR_black +
               NR_ivory , 
             data=all_class,surv=TRUE,x=TRUE, y=TRUE,time.inc=730)
surv <- Survival(f_cph)
nom_cph <- nomogram(f_cph, fun=list(function(x) surv(365, x),function(x) surv(730, x)),
                    fun.at = c(0.1,0.3,0.5,0.7,0.9),
                    funlabel=c("1-Year-Survival","2-Year-Survival"),
                    lp=F
)

jpeg('/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/nomogram_train.jpeg',
     width = 1200, height = 800, units = "px", pointsize = 12,quality = 100,bg = "#e5ecff",res=100)
plot(nom_cph, xfrac=.3,cex.axis=.7,cex.var=0.9,lmgp=.05,tcl=0.25) # delete R_white according to nomogram
dev.off()

cal2 <- calibrate(f_cph, cmethod='KM', method="boot", u=730, m=12)

jpeg('/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/calibration_curves_train.jpeg',
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





#test in SRP070710
read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/SRP070710_Symbol_FPKM_expr.txt",header = T,as.is = TRUE)->SRP070710_expr

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(SRA_Study == "SRP070710") %>%
  dplyr::filter(Survival_time != "NA")%>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(Run,Response,Survival_time,Survival_status,Age,Gender) ->metadata

list_sets[c("NR_black","NR_ivory")]->a
ssgava_score <- gsva(as.matrix(SRP070710_expr), a, min.sz=1, max.sz=999999, method="ssgsea",kcdf="Gaussian")

t(ssgava_score) %>%
  as.data.frame(stringsAsFactors=FALSE)->t_ssgava_score
cbind(rownames(t_ssgava_score),t_ssgava_score)->t_ssgava_score
rownames(t_ssgava_score)=NULL
colnames(t_ssgava_score)[1]="Run"

Combined_data=merge(metadata[,c(1,2,3,4)],t_ssgava_score)

# write.table(Combined_data,
#             "/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/SRP070710_selected.txt",
#             col.names = TRUE,row.names = FALSE,quote=FALSE,sep = "\t")

Combined_data->all_class

all_class%>%
  dplyr::filter(Response != "NE")->all_class

for (i in 1:nrow(all_class)) {
  if(all_class$Survival_status[i]=="Dead"){
    all_class$Survival_status[i]=2
  }else {
    all_class$Survival_status[i]=1
  }
}
all_class$Survival_status=as.numeric(all_class$Survival_status)
all_class$Survival_time=as.numeric(all_class$Survival_time)

ddist <- datadist(all_class)
options(datadist='ddist')

#cph(down order)
f_cph <- cph(Surv(Survival_time,Survival_status) ~
               Response +
               NR_black +
               NR_ivory, 
             data=all_class,surv=TRUE)
1-rcorrcens(Surv(Survival_time,Survival_status) ~ predict(f_cph), data =  all_class)#[[1]]

# Somers' Rank Correlation for Censored Data    Response variable:Surv(Survival_time, Survival_status)
# 
#                    C   Dxy  aDxy    SD     Z P   n
# predict(f_cph) 0.894 1.789 0.211 0.926 -9.68 1 -25

f_cph <- cph(Surv(Survival_time,Survival_status) ~
               Response +
               NR_black +
               NR_ivory , 
             data=all_class,surv=TRUE,x=TRUE, y=TRUE,time.inc=730)
# surv <- Survival(f_cph)
# nom_cph <- nomogram(f_cph, fun=list(function(x) surv(365, x),
#                                     function(x) surv(1095, x)),
#                     fun.at = c(0.1,0.3,0.5,0.7,0.9),
#                     funlabel=c("1-Year-Survival","3-Year-Survival"),
#                     lp=F
# )
# plot(nom_cph, xfrac=.3,cex.axis=.7,cex.var=0.8,lmgp=.05,tcl=0.25) # delete R_white according to nomogram
cal2 <- calibrate(f_cph, cmethod='KM', method="boot", u=730, m=8)

jpeg('/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/calibration_curves_test.jpeg',
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

write.table(a$NR_black,"/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/NR_black_genes.txt",
            row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
write.table(a$NR_ivory,"/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/NR_ivory_genes.txt",
            row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")

#enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)


read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/NR_black_genes.txt",sep="\t",
           header = F,as.is = TRUE) ->NR_black_genes
colnames(NR_black_genes)="Symbol"
read.table("/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/NR_ivory_genes.txt",sep="\t",
           header = F,as.is = TRUE) ->NR_ivory_genes
colnames(NR_ivory_genes)="Symbol"
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship

merge(relationship,NR_black_genes)->NR_black_genes
enrichGO(gene = NR_black_genes$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "fdr",pvalueCutoff = 1,readable = TRUE)->GO_enrich_black#
#DOSE::dotplot(GO_enrich_black, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")# p.adj>0.13
enrichKEGG(gene=NR_black_genes$GeneID,organism="human",pvalueCutoff=1,pAdjustMethod = "BH") ->KEGG_enrich_black# 0
enrichPathway(gene=NR_black_genes$GeneID,pvalueCutoff=1, readable=T)->Reactome_enrich_black# p.adj>0.12


merge(relationship,NR_ivory_genes)->NR_ivory_genes
enrichGO(gene = NR_ivory_genes$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "fdr",pvalueCutoff = 0.05,readable = TRUE)->GO_enrich_ivory#42
enrichKEGG(gene=NR_ivory_genes$GeneID,organism="human",pvalueCutoff=1,pAdjustMethod = "BH") ->KEGG_enrich_ivory# Ribosome p.adj=0.006224198 
enrichPathway(gene=NR_ivory_genes$GeneID,pvalueCutoff=0.05, readable=T)->Reactome_enrich_ivory

write.table(as.data.frame(GO_enrich_ivory),"/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/GO_enrich_ivory.txt",
            quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(as.data.frame(Reactome_enrich_ivory),"/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/Reactome_enrich_ivory.txt",
            quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)


union(NR_black_genes$GeneID,NR_ivory_genes$GeneID)->all_genes
enrichGO(gene = all_genes,OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "fdr",pvalueCutoff = 0.05,readable = TRUE)->GO_enrich_All#
enrichKEGG(gene=all_genes,organism="human",pvalueCutoff=1,pAdjustMethod = "BH") ->KEGG_enrich_All## p.adj>0.14
enrichPathway(gene=all_genes,pvalueCutoff=0.05, readable=T)->Reactome_enrich_All

write.table(as.data.frame(GO_enrich_All),"/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/GO_enrich_All.txt",
            quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
write.table(as.data.frame(Reactome_enrich_All),"/data/liull/immune-checkpoint-blockade/coexpress_modules/test_WGCNA_R/SRP094781/nomogram/Reactome_enrich_All.txt",
            quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
