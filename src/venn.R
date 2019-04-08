library (VennDiagram)

read.table("/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/up.txt",sep="\t",header = T,as.is = TRUE)->gastric_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/down.txt",sep="\t",header = T,as.is = TRUE)->gastric_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_up.txt",header = T,as.is = TRUE)->melanoma_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_down.txt",header = T,as.is = TRUE)->melanoma_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_up.txt",header = T,as.is = TRUE)->melanoma_CTLA4_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_down.txt",header = T,as.is = TRUE)->melanoma_CTLA4_down


venn.diagram(list(gastric_PD1=gastric_PD1_up$EnsemblId, melanoma_PD1=melanoma_PD1_up$EnsemblId,melanoma_CTLA4=melanoma_CTLA4_up$EnsemblId),filename=NULL,fill = c("cornflowerblue", "green", "yellow"))->venn_up
grid.draw(venn_up)

ggsave(
  filename = 'venn_up.pdf',
  plot = venn_up,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/different_expression/venn',
  width = 12,
  height = 8
)


venn.diagram(list(gastric_PD1=gastric_PD1_down$EnsemblId, melanoma_PD1=melanoma_PD1_down$EnsemblId,melanoma_CTLA4=melanoma_CTLA4_down$EnsemblId),filename=NULL,fill = c("cornflowerblue", "green", "yellow"))->venn_down
grid.draw(venn_down)

ggsave(
  filename = 'venn_down.pdf',
  plot = venn_down,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/different_expression/venn',
  width = 12,
  height = 8
)


# intersect(gastric_PD1_down$EnsemblId,melanoma_PD1_down$EnsemblId)
# #gastric_PD1&melanoma_PD1:MMP2,PIR,CYBRD1,COL6A1,COL3A1
# intersect(gastric_PD1_down$EnsemblId,melanoma_CTLA4_down$EnsemblId)
# #gastric_PD1&melanoma_CTLA4:NID1,ATP2B1,FSTL1,STXBP5,NDRG1
# intersect(melanoma_PD1_down$rowname,melanoma_CTLA4_down$rowname)
# #melanoma_PD1&melanoma_CTLA4:UBE2E3,MRC2 

intersect(melanoma_PD1_down$rowname,melanoma_CTLA4_down$rowname)
#CLEC18B  ,GDNF  ,
intersect(melanoma_PD1_up$rowname,melanoma_CTLA4_up$rowname)
#CLEC17A  ,DAZ1  ,DAZ4  ,CR2  ,DAZ2


