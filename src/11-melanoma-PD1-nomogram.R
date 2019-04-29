read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/R_skyblue3_class.txt",
           header = T,as.is = TRUE)->R_skyblue3_class
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/NR_midnightblue_class.txt",
           header = T,as.is = TRUE)->NR_midnightblue_class
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/NR_orangered4_class.txt",
           header = T,as.is = TRUE)->NR_orangered4_class
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/NR_white_class.txt",
           header = T,as.is = TRUE)->NR_white_class
# read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/8_gene_class.txt",
#            header = T,as.is = TRUE)->eight_gene_class
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

#psm(up order)-----------------------------------------------------------------------------------------------------------------------------------------
f_psm <- psm(Surv(Survival_time,Survival_status) ~ 
               Response +
               R_skyblue3_class +
               NR_midnightblue_class + 
               NR_orangered4_class + 
               NR_white_class, 
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



#cph(down order)---------------------------------
f_cph <- cph(Surv(Survival_time,Survival_status) ~
               Response +
               R_skyblue3_class +
               NR_midnightblue_class +
               NR_orangered4_class +
               NR_white_class,
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

#predicted~actual------------------
rcorrcens(Surv(Survival_time,Survival_status) ~ predict(f_psm), data =  all_class)#0.811
rcorrcens(Surv(Survival_time,Survival_status) ~ predict(f_cph), data =  all_class)#0.187
###(psm)---------------------------

# f2 <- psm(Surv(Survival_time,Survival_status) ~ 
#             R_skyblue3_class + 
#             NR_midnightblue_class + 
#             NR_orangered4_class + 
#             NR_white_class, 
#           data =  all_class, x=T, y=T, dist='lognormal',time.inc=365) ===f_psm

cal1 <- calibrate(f_psm, cmethod='KM', u=365, m=16,B=70)

jpeg('/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/psm_validation.jpeg',
     width = 800, height = 800, units = "px", pointsize = 12,quality = 100,bg = "#e5ecff",res=100)
par(mar=c(5,5,3,3),cex = 1.0)
plot(cal1,xlim=c(0,1),ylim=c(0,1),subtitles=FALSE, 
     cex.subtitles=.75, riskdist=FALSE)
abline(0,1,lty =3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))
dev.off()

###(cph)---------------------------

# f <- cph(Surv(Survival_time,Survival_status) ~ 
#            R_skyblue3_class + 
#            NR_midnightblue_class + 
#            NR_orangered4_class + 
#            NR_white_class, 
#          x=TRUE, y=TRUE,surv = TRUE,time.inc=365,data=all_class)  ===f_cph

cal<-calibrate(f_cph,u=365,cmethod='KM',m=16,B=70)

jpeg('/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/WGCNA/cph_validation.jpeg',
     width = 800, height = 800, units = "px", pointsize = 12,quality = 100,bg = "#e5ecff",res=100)
par(mar=c(5,5,3,3),cex = 1.0)
plot(cal,xlim = c(0,1),ylim= c(0,1),subtitles=FALSE, 
     cex.subtitles=.75, riskdist=FALSE)
abline(0,1,lty =3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))
dev.off()



