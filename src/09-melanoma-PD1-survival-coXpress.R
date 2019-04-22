library(coXpress)
library(magrittr)


# Step 1: reading your data
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/pre_PD1_filtered_symbol_expr.txt",header = T,as.is = TRUE) ->pre_PD1_expr

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response,Biopsy_Time,Survival_time,Survival_status,Age,Gender) ->metadata
dplyr::filter(metadata,Response %in% c("CR","PR","R"))-> response#26
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response#59
dplyr::select(pre_PD1_expr,response$Run,non_response$Run)->ordered_PD1_expr


#coXpress based on Response
hc.gene  <- cluster.gene(ordered_PD1_expr[,1:26],s="pearson",m="average")
g <- cutree(hc.gene, h=0.4)
cox <- coXpress(ordered_PD1_expr, g, 1:26, 27:85,times=1000)
plot.compare.group(ordered_PD1_expr,g,21,1:26,27:85,
                   scale.center=TRUE,scale.scale=TRUE,
                   ylim=c(-5,5))
inspect.group(ordered_PD1_expr, g, 21, 1:26, 27:85)

cox[cox$N>8 & cox$pr.g1<=0.05 & cox$pr.g2>=0.05, ]->selected_group_R
cox[cox$N>8 & cox$pr.g1>0.05 & cox$pr.g2<0.05, ]->selected_group_NR
rbind(selected_group_R,selected_group_NR)->selected_groups

write.table(g,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/g.txt",quote = FALSE,col.names = FALSE,row.names = TRUE)
write.table(cox,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/all_groups.txt",quote = FALSE,col.names = TRUE)
write.table(selected_groups,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/selected_groups.txt",quote = FALSE,col.names = TRUE)





#coXpress R 1000 times

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/selected_groups.txt",header = T,as.is = TRUE)->selected_groups

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/g.txt",header = F,as.is = TRUE) ->g_dataframe
g=g_dataframe$V2
names(g)=g_dataframe$V1

list_sets=list()
module_names=character()
for (i in 1:nrow(selected_groups)) {
  
  module=selected_groups$group[i]
  inspect.group(ordered_PD1_expr,g,selected_groups$group[i],1:26, 27:85)->sets
  list_sets[i]=list(union(sets$GeneA,sets$GeneB))
  module_names[i]=paste("R_module",module,sep="_")
  
}
names(list_sets)=module_names



#ssGSEA--------------------

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Survival_time != "NA")%>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%###
  dplyr::select(Run,Response,Survival_time,Survival_status,Age,Gender) ->metadata

ssgava_score <- gsva(as.matrix(pre_PD1_expr), list_sets, min.sz=1, max.sz=999999, method="ssgsea",kcdf="Gaussian", abs.ranking=FALSE, verbose=TRUE)

t(ssgava_score) %>%
  as.data.frame(stringsAsFactors=FALSE)->t_ssgava_score
cbind(rownames(t_ssgava_score),t_ssgava_score)->t_ssgava_score
rownames(t_ssgava_score)=NULL
colnames(t_ssgava_score)[1]="Run"

Combined_data=merge(metadata[,c(1,3,4,5,6)],t_ssgava_score)

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
colnames(filter_res)[1]="gene_module"
filter_res %>%
  dplyr::filter(p.value<=0.05)%>%
  dplyr::select(gene_module)->selected_modules
# gene_module
# 1   R_module_239


Combined_data %>%
  dplyr::select(Run,Survival_time,Survival_status,as.character(selected_modules[1,1]))%>%
  dplyr::mutate(Class=rep("class",nrow(Combined_data)))-> Combined_module_1
cutoff=mean(Combined_module_1[,4])

for (i in 1:nrow(Combined_module_1)) {
  if(Combined_module_1[i,4]>=cutoff){
    Combined_module_1$Class[i]="high"
  }else {
    Combined_module_1$Class[i]="low"
  }
}

fit <- survfit(Surv(Survival_time, Survival_status) ~ Class, data = Combined_module_1)

ggsurvplot(fit, data = Combined_module_1, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")#0.51


write.table(list_sets$R_module_239,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/R_module_239.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
write.table(list_sets$R_module_1448,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/R_module_1448.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
write.table(list_sets$R_module_40,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/R_module_40.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
write.table(list_sets$R_module_51,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/R_module_51.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
write.table(list_sets$R_module_1271,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/R_module_1271.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
write.table(list_sets$R_module_113,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/R_module_113.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")
write.table(list_sets$R_module_294,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/R_module_294.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")

