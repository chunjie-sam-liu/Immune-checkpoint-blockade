library("magrittr")
library("ggplot2")
library("ggpubr")
library("survminer")
library("survival")

fileName<-dir("E:/OneDrive/GitHub/tumor immunology/survplot/50%")

for (i in 1:24){
  yourfilename=paste(fileName[i],".pdf",sep="")
  pdf(file=yourfilename)
  Surplot <- read.csv(file = fileName[i], header = TRUE , check.names = FALSE)
  fit <- survfit(Surv(suppressWarnings(as.numeric(as.character(Survival_time))), as.numeric(Survival_status)) ~ type, data = Surplot)
  print(i)
  print(ggsurvplot(fit, data = Surplot,pval = TRUE, xlab = "Time in days"))

  dev.off()
}