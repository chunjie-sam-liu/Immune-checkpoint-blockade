library("magrittr")
library("ggplot2")
library("ggpubr")
library("survminer")
library("survival")

Surplot <- read.csv(file = "24_CD8_T.csv" , header = TRUE , check.names = FALSE)
fit <- survfit(Surv(suppressWarnings(as.numeric(as.character(Survival_time))), as.numeric(Survival_status)) ~ type, data = Surplot)

ggsurvplot(fit, data = Surplot,pval = TRUE, xlab = "Time in days")
