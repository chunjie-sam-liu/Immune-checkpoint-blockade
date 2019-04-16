library (VennDiagram)

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_ENSG.txt",header = T,as.is = TRUE)->gastric_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_ENSG.txt",header = T,as.is = TRUE)->gastric_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_up_ENSG.txt",header = T,as.is = TRUE)->melanoma_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_down_ENSG.txt",header = T,as.is = TRUE)->melanoma_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_up_ENSG.txt",header = T,as.is = TRUE)->melanoma_CTLA4_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_down_ENSG.txt",header = T,as.is = TRUE)->melanoma_CTLA4_down


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


intersect(gastric_PD1_up$Ensembl_ID,melanoma_PD1_up$rowname)
#IDO1, JAKMIP1,TRAT1,LRG1,LCK ,NELL2,PDCD1,HLA-DQA1 ,HLA-DOA ,HCP5,UBD ,LGALS17A,OR2I1P
#"NONHSAG213716"  "NONHSAG223513"  "NONHSAG223686" "NONHSAG224412" "NONHSAG268497" "NONHSAG268538" "NONHSAG270697""NONHSAG280970"
intersect(gastric_PD1_down$Symbol,melanoma_PD1_down$Symbol)
#VCAN,VEGFD,GFRA3,MMP2,NTF3,ROS1,CORO2B,OLFM1,RGS22,NSG1,TM6SF2,MUC13,ANO2,HPSE2,HIF3A,ADAM33,MAGED4B,MS4A8,ESYT3,CCDC8,SCUBE3,RCOR2,RSPO4,ZNF662,MAGED4
#"NONHSAG229342"   "NONHSAG229343"   "NONHSAG229344"  "NONHSAG260990"   "NONHSAG265471"   "NONHSAG284812" 

