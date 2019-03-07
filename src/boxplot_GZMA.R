library(magrittr)

#1
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_PD1_removed_batch_expression.txt",header = T,as.is = TRUE) ->melanoma_PD1
cbind(rownames(melanoma_PD1),melanoma_PD1)->melanoma_PD1
colnames(melanoma_PD1)[1]="ensembl_ID"

dplyr::filter(melanoma_PD1,ensembl_ID=="ENSG00000145649")->GZMA_1

readxl::read_excel("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata
dplyr::filter(metadata,Response %in% c("CR","PR","PRCR","R")) -> response
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response

group_1=c(rep("CR/PR",nrow(response)),rep("SD/PD",nrow(non_response)))
dplyr::select(GZMA_1,response$Run,non_response$Run)%>%
  rbind(group_1)%>%
  t()%>%
  as.data.frame()->ordered_GZMA_1
colnames(ordered_GZMA_1)=c("GZMA","group")
ordered_GZMA_1$GZMA=as.numeric(as.character(ordered_GZMA_1$GZMA))
ordered_GZMA_1$GZMA[which(ordered_GZMA_1$GZMA[]<0)]=0.001
ordered_GZMA_1$GZMA=log2(ordered_GZMA_1$GZMA)
#ggboxplot(ordered_GZMA_1, x="group", y="GZMA", color = "group")+stat_compare_means()

#2
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_CTLA4_removed_batch_expression.txt",header = T,as.is = TRUE) ->melanoma_CTLA4
cbind(rownames(melanoma_CTLA4),melanoma_CTLA4)->melanoma_CTLA4
colnames(melanoma_CTLA4)[1]="ensembl_ID"

dplyr::filter(melanoma_CTLA4,ensembl_ID=="ENSG00000145649")->GZMA_2

readxl::read_excel("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") -> SRA
readxl::read_excel("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP") -> dbGAP  
rbind(SRA,dbGAP) %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata
dplyr::filter(metadata,Response %in% c("CR","PR","PRCR","R")) -> response
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response

group_2=c(rep("CR/PR",nrow(response)),rep("SD/PD",nrow(non_response)))
dplyr::select(GZMA_2,response$Run,non_response$Run)%>%
  rbind(group_2)%>%
  t()%>%
  as.data.frame()->ordered_GZMA_2
colnames(ordered_GZMA_2)=c("GZMA","group")
ordered_GZMA_2$GZMA=as.numeric(as.character(ordered_GZMA_2$GZMA))
#ordered_GZMA_2$GZMA[which(ordered_GZMA_2$GZMA[]<0)]=0.001#no less than 0
ordered_GZMA_2$GZMA=log2(ordered_GZMA_2$GZMA)
#ggboxplot(ordered_GZMA_2, x="group", y="GZMA", color = "group")+stat_compare_means()

#3
readxl::read_excel("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="gastric cancer") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

rbind(dplyr::filter(metadata,Response=="CR"),dplyr::filter(metadata,Response=="PR")) ->response
rbind(dplyr::filter(metadata,Response=="SD"),dplyr::filter(metadata,Response=="PD")) ->non_response


read.table("/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/gastric_cancer_PD1_DEG.txt",sep="\t",header = T,as.is = TRUE) %>%
  dplyr::select(ensembl_ID,response$Run,non_response$Run)%>%
  dplyr::filter(ensembl_ID=="ENSG00000145649")->GZMA_3

group_3=c(rep("CR/PR",nrow(response)),rep("SD/PD",nrow(non_response)))
GZMA_3[,-1]%>%
  rbind(group_3)%>%
  t()%>%
  as.data.frame()->ordered_GZMA_3
colnames(ordered_GZMA_3)=c("GZMA","group")
ordered_GZMA_3$GZMA=log2(as.numeric(as.character(ordered_GZMA_3$GZMA)))
#ggboxplot(ordered_GZMA_3, x="group", y="GZMA", color = "group")+stat_compare_means(method = "wilcox.test")

#plot
cbind(rep("melanoma_PD1",nrow(ordered_GZMA_1)),ordered_GZMA_1)->a
cbind(rep("melanoma_CTLA4",nrow(ordered_GZMA_2)),ordered_GZMA_2)->b
cbind(rep("gastric_cancer_PD1",nrow(ordered_GZMA_3)),ordered_GZMA_3)->c
colnames(a)[1]="Cancer"
colnames(b)[1]="Cancer"
colnames(c)[1]="Cancer"
rbind(a,b,c)->all

anno1 <- wilcox.test(ordered_GZMA_1[ordered_GZMA_1$group=="CR/PR",1],ordered_GZMA_1[ordered_GZMA_1$group=="SD/PD",1])$p.value
anno2 <- wilcox.test(ordered_GZMA_2[ordered_GZMA_2$group=="CR/PR",1],ordered_GZMA_2[ordered_GZMA_2$group=="SD/PD",1])$p.value
anno3 <- wilcox.test(ordered_GZMA_3[ordered_GZMA_3$group=="CR/PR",1],ordered_GZMA_3[ordered_GZMA_3$group=="SD/PD",1])$p.value
ggplot(all, aes(x=Cancer, y=GZMA, fill = group))+
  geom_boxplot(position="dodge")+
  theme(panel.grid =element_blank())+
  geom_signif(annotation=formatC(anno1, digits=3),y_position=7.5, xmin=0.8, xmax=1.19, margin_top= 0.05)+
  geom_signif(annotation=formatC(anno2, digits=3),y_position=7.5, xmin=1.8, xmax=2.19)+
  geom_signif(annotation=formatC(anno3, digits=3),y_position=7.5, xmin=2.8, xmax=3.19) -> p

#p <- ggboxplot(all, x="Cancer", y="GZMA", color = "group", palette = "jco") + 
#  stat_compare_means(aes(group=group,label=..p.format..),label.y = c(9, 9, 9))
ggsave(
  filename = 'boxplot_GZMA_3.pdf',
  plot = p,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/different_expression',
  width = 6,
  height = 6.8
)
