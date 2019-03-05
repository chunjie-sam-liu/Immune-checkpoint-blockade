read.table("/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/up.txt",sep="\t",header = T,as.is = TRUE)->gastric_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/down.txt",sep="\t",header = T,as.is = TRUE)->gastric_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/up.txt",sep="\t",header = T,as.is = TRUE)->melanoma_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/down.txt",sep="\t",header = T,as.is = TRUE)->melanoma_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/up.txt",sep="\t",header = T,as.is = TRUE)->melanoma_CTLA4_up
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/down.txt",sep="\t",header = T,as.is = TRUE)->melanoma_CTLA4_down

gastric_PD1=c(gastric_PD1_up$EnsemblId,gastric_PD1_down$EnsemblId)
melanoma_PD1=c(melanoma_PD1_up$EnsemblId,melanoma_PD1_down$EnsemblId)
melanoma_CTLA4=c(melanoma_CTLA4_up$EnsemblId,melanoma_CTLA4_down$EnsemblId)
