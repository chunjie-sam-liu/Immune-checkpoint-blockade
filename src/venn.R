library (VennDiagram)

read.table("/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/up.txt",sep="\t",header = T,as.is = TRUE)->gastric_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/down.txt",sep="\t",header = T,as.is = TRUE)->gastric_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/up.txt",sep="\t",header = T,as.is = TRUE)->melanoma_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/down.txt",sep="\t",header = T,as.is = TRUE)->melanoma_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/up.txt",sep="\t",header = T,as.is = TRUE)->melanoma_CTLA4_up
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/down.txt",sep="\t",header = T,as.is = TRUE)->melanoma_CTLA4_down


venn.diagram(list(gastric_PD1=gastric_PD1_up$EnsemblId, melanoma_PD1=melanoma_PD1_up$EnsemblId,melanoma_CTLA4=melanoma_CTLA4_up$EnsemblId),filename=NULL,fill = c("cornflowerblue", "green", "yellow"))->venn_up
grid.draw(venn_up)
#CD74,GZMA,CCL5

venn.diagram(list(gastric_PD1=gastric_PD1_down$EnsemblId, melanoma_PD1=melanoma_PD1_down$EnsemblId,melanoma_CTLA4=melanoma_CTLA4_down$EnsemblId),filename=NULL,fill = c("cornflowerblue", "green", "yellow"))->venn_down
grid.draw(venn_down)
#gastric_PD1&melanoma_PD1:TNC,MMP2,MFAP2,MYO10    SERPINE2,PPFIBP1,PIR,BUD31,DNAJC13,PHF20L1,PIP4P2,GADD45GIP1
#gastric_PD1&melanoma_CTLA4:ELN,NID1,GJA1   ATP2B1,GFPT1,NR3C1,CNOT3,PPFIBP2,NDRG1,FSTL1,STXBP5,TRIM74
#melanoma_PD1&melanoma_CTLA4:UBE2E3,AP3M2,GLOD4,MRPS6 

