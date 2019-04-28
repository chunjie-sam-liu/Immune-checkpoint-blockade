read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/R_skyblue3_class.txt",
           header = T,as.is = TRUE)->R_skyblue3_class
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/NR_midnightblue_class.txt",
           header = T,as.is = TRUE)->NR_midnightblue_class
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/NR_orangered4_class.txt",
           header = T,as.is = TRUE)->NR_orangered4_class
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/NR_white_class.txt",
           header = T,as.is = TRUE)->NR_white_class
merge(R_skyblue3_class,NR_midnightblue_class)%>%
  merge(NR_orangered4_class)%>%
  merge(NR_white_class)->all_class

for (i in 1:nrow(all_class)) {
  if(all_class$Survival_status[i]=="Dead"){
    all_class$Survival_status[i]=2
  }else {
    all_class$Survival_status[i]=1
  }
}

all_class%>%
  dplyr::filter(Response != "NE")->all_class

all_class$Survival_status=as.numeric(all_class$Survival_status)
all_class$Survival_time=as.numeric(all_class$Survival_time)

ddist <- datadist(all_class)
options(datadist='ddist')

#psm-----------------------------------------------------------------------------------------------------------------------------------------
f <- psm(Surv(Survival_time,Survival_status) ~ Response + R_skyblue3_class + NR_midnightblue_class + NR_orangered4_class + NR_white_class, 
         data=all_class, dist='lognormal')
med  <- Quantile(f)
surv <- Survival(f)  # This would also work if f was from cph
plot(nomogram(f, fun=function(x) med(lp=x), funlabel="Median Survival Time"))

nom <- nomogram(f, fun=list(function(x) surv(1095, x),
                            function(x) surv(1825, x)),
                funlabel=c("1-Year Survival Probability", 
                           "3-Year Survival Probability"))
plot(nom, xfrac=.7)



# cox_f <- cph(Surv(Survival_time,Survival_status) ~ Response + R_skyblue3_class + NR_midnightblue_class + NR_orangered4_class + NR_white_class, 
#          data=all_class,surv=TRUE)
# surv_cox <- Survival(cox_f)
# nom_cox <- nomogram(cox_f, fun=list(function(x) surv_cox(1095, x),
#                                     function(x) surv_cox(1825, x)),
#                     funlabel=c("1-Year Survival Probability", 
#                                "3-Year Survival Probability"),
#                     lp=F)
#function(x) surv.cox(1825, x)           ,"5-Year-Survival"

#cox----------------------------------
mod.cox=cph(Surv(Survival_time,Survival_status)~R_skyblue3_class+NR_midnightblue_class+NR_orangered4_class+NR_white_class,all_class,surv=TRUE)

surv.cox <- Survival(mod.cox)
nom.cox <- nomogram(mod.cox, fun=list(function(x) surv.cox(365, x),
                                      function(x) surv.cox(1095, x)),
                    funlabel=c("1-Year-Survival","3-Year-Survival"),
                    lp=F)
plot(nom.cox)
