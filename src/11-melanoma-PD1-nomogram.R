read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/R_skyblue3_class.txt",
           header = T,as.is = TRUE)->R_skyblue3_class
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/NR_midnightblue_class.txt",
           header = T,as.is = TRUE)->NR_midnightblue_class
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/NR_orangered4_class.txt",
           header = T,as.is = TRUE)->NR_orangered4_class
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/NR_white_class.txt",
           header = T,as.is = TRUE)->NR_white_class
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/8_gene_class.txt",
           header = T,as.is = TRUE)->eight_gene_class
read.table("/data/liull/immune-checkpoint-blockade/TIL/melanoma_PD1/TIL_classify.txt",
           header = T,as.is = TRUE) %>% dplyr::select(Run,Tr1)-> Tr1_class

merge(R_skyblue3_class,NR_midnightblue_class)%>%
  merge(NR_orangered4_class)%>%
  merge(NR_white_class) %>%
  merge(Tr1_class)%>%
  merge(eight_gene_class)->all_class

all_class%>%
  dplyr::filter(Response != "NE")->all_class


#Survival_time,Survival_status ~ modules---------------------------------------------------
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

#psm(up order)
f_psm <- psm(Surv(Survival_time,Survival_status) ~ 
               Response +
               R_skyblue3_class +
               NR_midnightblue_class + 
               NR_orangered4_class + 
               NR_white_class + Tr1 + eight_gene_class, 
             data=all_class, x=T, y=T, dist='lognormal',time.inc=365)
med  <- Quantile(f_psm)
surv <- Survival(f_psm)  # This would also work if f was from cph
plot(nomogram(f_psm, fun=function(x) med(lp=x), funlabel="Median Survival Time"))

nom_psm <- nomogram(f_psm, fun=list(function(x) surv(365, x),
                            function(x) surv(1095, x)),
                    fun.at = c(0.1,0.2,0.4,0.6,0.8,0.9,1),
                    funlabel=c("1-Year Survival Probability", 
                           "3-Year Survival Probability"))

jpeg('/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/psm_nomogram.jpeg',
     width = 1000, height = 500, units = "px", pointsize = 12,quality = 100,bg = "#e5ecff",res=100)
plot(nom_psm, xfrac=.4,cex.axis=.7)
dev.off()



#cph(down order)
f_cph <- cph(Surv(Survival_time,Survival_status) ~
               Response +
               R_skyblue3_class +
               NR_midnightblue_class +
               NR_orangered4_class +
               NR_white_class +
               Tr1,
             data=all_class,surv=TRUE,x=TRUE, y=TRUE,time.inc=365)

surv <- Survival(f_cph)
nom_cph <- nomogram(f_cph, fun=list(function(x) surv(365, x),
                                    function(x) surv(1095, x)),
                    fun.at = c(0.1,0.3,0.5,0.7,0.9),
                    funlabel=c("1-Year-Survival","3-Year-Survival"),
                    lp=F
                    )
jpeg('/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/cph_nomogram.jpeg',
     width = 1000, height = 500, units = "px", pointsize = 12,quality = 100,bg = "#e5ecff",res=100)
plot(nom_cph, xfrac=.3,cex.axis=.7,cex.var=0.8,lmgp=.05,tcl=0.25)
dev.off()

#actual ~ predict
rcorrcens(Surv(Survival_time,Survival_status) ~ predict(f_psm), data =  all_class)#0.811
rcorrcens(Surv(Survival_time,Survival_status) ~ predict(f_cph), data =  all_class)#0.187
# rcorrcens(Surv(Survival_time,Survival_status) ~ predict(f_psm), data =  all_class)#NR~R 0.789
# rcorrcens(Surv(Survival_time,Survival_status) ~ predict(f_cph), data =  all_class)#NR~R 0.211

###(psm)

cal1 <- calibrate(f_psm, cmethod='KM', u=365, m=16,B=70)

jpeg('/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/psm_validation.jpeg',
     width = 800, height = 800, units = "px", pointsize = 12,quality = 100,bg = "#e5ecff",res=100)
par(mar=c(5,5,3,3),cex = 1.0)
plot(cal1,xlim=c(0,1),ylim=c(0,1),subtitles=FALSE, 
     cex.subtitles=.75, riskdist=FALSE)
abline(0,1,lty =3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))
dev.off()

###(cph)

cal<-calibrate(f_cph,u=365,cmethod='KM',m=16,B=70)

jpeg('/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/cph_validation.jpeg',
     width = 800, height = 800, units = "px", pointsize = 12,quality = 100,bg = "#e5ecff",res=100)
par(mar=c(5,5,3,3),cex = 1.0)
plot(cal,xlim = c(0,1),ylim= c(0,1),subtitles=FALSE, 
     cex.subtitles=.75, riskdist=FALSE)
abline(0,1,lty =3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))
dev.off()

#response ~ modules---------------------------------------------------------------------
all_class$Response %>% 
  gsub("^CR$",1,.)%>%
  gsub("^PR$",1,.)%>%
  gsub("^SD$",0,.)%>%
  gsub("^PD$",0,.)->all_class$Response


all_class$Response=as.numeric(all_class$Response)

ddist <- datadist(all_class)
options(datadist='ddist')

fit <- lrm(Response ~ NR_midnightblue_class +  NR_white_class + Tr1 + eight_gene_class,
           data=all_class,x=TRUE,y=TRUE)


nom <- nomogram(fit, fun=plogis,  # or fun=plogis
                funlabel="Possibility of Response")
jpeg('/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/R_NR_nomogram.jpeg',
     width = 800, height = 800, units = "px", pointsize = 12,quality = 100,bg = "#e5ecff",res=100)
plot(nom, xfrac=.45)
dev.off()

rcorrcens(Response ~ predict(fit), data =  all_class)
# Somers' Rank Correlation for Censored Data    Response variable:Response
# 
#                  C   Dxy  aDxy    SD   Z P  n
# predict(fit) 0.789 0.578 0.578 0.109 5.3 0 75
