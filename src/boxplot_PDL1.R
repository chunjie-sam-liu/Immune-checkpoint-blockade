library(magrittr)

read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_PD1_removed_batch_expression.txt",header = T,as.is = TRUE) ->melanoma_PD1
cbind(rownames(melanoma_PD1),melanoma_PD1)->melanoma_PD1
colnames(melanoma_PD1)[1]="ensembl_ID"

dplyr::filter(melanoma_PD1,ensembl_ID=="ENSG00000120217")->PDL1


readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

#select response and non_response's sample id and project id
dplyr::filter(metadata,Response %in% c("CR","PR","PRCR","R")) -> response
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response

group=c(rep("CR/PR",nrow(response)),rep("SD/PD",nrow(non_response)))
dplyr::select(PDL1,response$Run,non_response$Run) %>%
  rbind(group)%>%
  t()%>%
  as.data.frame()->ordered_PDL1
colnames(ordered_PDL1)=c("PDL1","group")
ordered_PDL1$PDL1=as.numeric(ordered_PDL1$PDL1)

ggplot(ordered_PDL1,aes(x=group,y=PDL1,fill=group))+geom_boxplot()+labs(title="melanoma_PD1")->boxplot_melanoma_PDL1
ggsave(
  filename = 'boxplot_melanoma_PD1_PDL1.pdf',
  plot = boxplot_melanoma_PDL1,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/different_expression',
  width = 6,
  height = 6.8
)
#2-----------------------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_CTLA4_removed_batch_expression.txt",header = T,as.is = TRUE) ->melanoma_CTLA4
cbind(rownames(melanoma_CTLA4),melanoma_CTLA4)->melanoma_CTLA4
colnames(melanoma_CTLA4)[1]="ensembl_ID"

dplyr::filter(melanoma_CTLA4,ensembl_ID=="ENSG00000120217")->PDL1


readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") -> SRA
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP") -> dbGAP  
  rbind(SRA,dbGAP) %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

#select response and non_response's sample id and project id
dplyr::filter(metadata,Response %in% c("CR","PR","PRCR","R")) -> response
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response

group=c(rep("CR/PR",nrow(response)),rep("SD/PD",nrow(non_response)))
dplyr::select(PDL1,response$Run,non_response$Run) %>%
  rbind(group)%>%
  t()%>%
  as.data.frame()->ordered_PDL1
colnames(ordered_PDL1)=c("PDL1","group")
ordered_PDL1$PDL1=as.numeric(ordered_PDL1$PDL1)

ggplot(ordered_PDL1,aes(x=group,y=PDL1,fill=group))+geom_boxplot()+labs(title="melanoma_CTLA4")->boxplot_melanoma_CTLA4
ggsave(
  filename = 'boxplot_melanoma_CTLA4_PDL1.pdf',
  plot = boxplot_melanoma_CTLA4,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/different_expression',
  width = 6,
  height = 6.8
)

#3----------------------------------------------------------------

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="gastric cancer") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(SRA_Study,Run,Response) ->metadata

rbind(dplyr::filter(metadata,Response=="CR"),dplyr::filter(metadata,Response=="PR")) ->response
rbind(dplyr::filter(metadata,Response=="SD"),dplyr::filter(metadata,Response=="PD")) ->non_response


read.table("/data/liull/immune-checkpoint-blockade/expression/all_FPKM_expression_2.txt",sep="\t",header = T,as.is = TRUE) %>%
  dplyr::select(gene_id,response$Run,non_response$Run)%>%
  dplyr::filter(gene_id=="ENSG00000120217")->PDL1

group=c("group",rep("CR/PR",nrow(response)),rep("SD/PD",nrow(non_response)))
rbind(PDL1,group)->ordered_PDL1
rownames(ordered_PDL1)=c("PDL1","group")
ordered_PDL1=ordered_PDL1[,-1]
ordered_PDL1=as.data.frame(t(ordered_PDL1))
ordered_PDL1$PDL1=as.numeric(ordered_PDL1$PDL1)
ggplot(ordered_PDL1,aes(x=group,y=PDL1,fill=group))+geom_boxplot()+labs(title="gastric_cancer_PD1")->boxplot_gastric_cancer_PDL1
ggsave(
  filename = 'boxplot_gastric_cancer_PD1_PDL1.pdf',
  plot = boxplot_gastric_cancer_PDL1,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/different_expression',
  width = 6,
  height = 6.8
)
