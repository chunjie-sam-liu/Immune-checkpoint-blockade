library(coXpress)
library(magrittr)


# Step 1: reading your data
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/pre_PD1_filtered_symbol_expr.txt",header = T,as.is = TRUE) ->all_expr

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response,Biopsy_Time,Survival_time,Survival_status,Age,Gender) ->metadata
dplyr::filter(metadata,Response %in% c("CR","PR","R"))-> response
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response
dplyr::select(all_expr,response$Run,non_response$Run)->ordered_all_expr

#focus on response group---------------------------------------------------------------------------------------------------------------------------
# Step 2: Cluster data based on a subset of experiments.
hc.gene_R <- cluster.gene(ordered_all_expr[,1:nrow(response)], "pearson", "average")

# Step 3: cut the tree at a height of 0.4 
g_R <- cutree(hc.gene_R, h=0.4)

# Step 4: examine the difference between R and NR
cox_R <- coXpress(ordered_all_expr, g_R, 1:nrow(response), (nrow(response)+1):(nrow(response)+nrow(non_response)), times=1000)
write.table(cox_R,"all_groups_R.txt",quote = FALSE,col.names = TRUE)

# Step 5: the results are a data.frame, with one row for each group of genes.
cox_R[cox_R$pr.g1 <= 0.05 & cox_R$pr.g2 >= 0.05 & cox_R$N>=7,]->selected_group_R

# Step 6: examine groups of interest graphically ,look at group 239
plot.compare.group(ordered_all_expr,g,239,1:nrow(response), (nrow(response)+1):(nrow(response)+nrow(non_response)),
                   scale.center=TRUE,scale.scale=TRUE,
                   ylim=c(-5,5))
inspect.group(ordered_all_expr,g,239,1:nrow(response), (nrow(response)+1):(nrow(response)+nrow(non_response)))#group 239's paied members' correlation in the R and NR



#focus on non_response group---------------------------------------------------------------------------------------------------------
hc.gene_NR <- cluster.gene(ordered_all_expr[,(nrow(response)+1):(nrow(response)+nrow(non_response))], "pearson", "average")
g_NR <- cutree(hc.gene_NR, h=0.4)
cox_NR <- coXpress(ordered_all_expr, g_NR, (nrow(response)+1):(nrow(response)+nrow(non_response)),1:nrow(response),  times=1000)
write.table(cox_NR,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/all_groups_NR.txt",quote = FALSE,col.names = TRUE)

cox_NR[cox_NR$pr.g1 <= 0.05 & cox_NR$pr.g2 >= 0.05 & cox_NR$N>=7,]->selected_group_NR




#ssGSEA--------------------
list_sets=list()
group_names=character()
for (i in 1:nrow(selected_group)) {
  
  group=selected_group$group[i]
  inspect.group(ordered_all_expr,g,selected_group$group[i],1:26, 27:85)->sets
  list_sets[i]=list(union(sets$GeneA,sets$GeneB))
  group_names[i]=paste("group",group,sep="_")
  
}
names(list_sets)=group_names

ssgava_score <- gsva(as.matrix(ordered_all_expr), list_sets, min.sz=1, max.sz=999999, method="ssgsea",kcdf="Gaussian", abs.ranking=FALSE, verbose=TRUE)

t(ssgava_score) %>%
  as.data.frame(stringsAsFactors=FALSE)->t_ssgava_score
cbind(rownames(t_ssgava_score),t_ssgava_score)->t_ssgava_score
rownames(t_ssgava_score)=NULL
colnames(t_ssgava_score)[1]="Run"

Combined_data=merge(metadata[,c(2,5,6,7,8)],t_ssgava_score)

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
colnames(filter_res)[1]="gene_group"
filter_res %>%
  dplyr::filter(p.value<=0.05)%>%
  dplyr::select(gene_group)->groups
# > groups
# gene_group
# 1  group_239

Combined_data %>%
  dplyr::select(Run,Survival_time,Survival_status,as.character(groups[,1]))->Combined_data_2
cutoff=median(Combined_data_2$group_239)
Combined_data_2 %>%
  dplyr::mutate(Class=rep("class",nrow(Combined_data_2)))%>%
  dplyr::select(Run,Survival_time,Survival_status,group_239,Class)->Combined_data_3



for (i in 1:nrow(Combined_data_3)) {
  if(Combined_data_3$group_239[i]>=cutoff){
    Combined_data_3$Class[i]="high"
  }else {
    Combined_data_3$Class[i]="low"
  }
}

fit <- survfit(Surv(Survival_time, Survival_status) ~ Class, data = Combined_data_3)
ggsurvplot(fit, data = Combined_data_3, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")






