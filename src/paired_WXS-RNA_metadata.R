library(magrittr)

readxl::read_excel("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP")->dbGAP
readxl::read_excel("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA")->SRA
readxl::read_excel("/data/liucj/data/immune-checkpoint-blockade/paired_metadata.xlsx",col_names = TRUE)->paired_metadata
colnames(paired_metadata)[2]="WXS_Status"

#SRP011540 WXS_normal~WXS_tumor~RNA_Seq
dbGAP %>%
  dplyr::filter(SRA_Study=="SRP011540")%>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::select(Run,Biopsy_Time,patient_id) -> SRP011540
colnames(SRP011540)=c("RNA_seq_tumor","RNA_seq_Status","Patient_id")

paired_metadata %>%
  dplyr::filter(Project=="SRP011540")%>%
  merge(SRP011540,all=F) %>%
  dplyr::select(Cancer,Immune_checkpoint,Project,Patient_id,Response,Second_Response_standard,WXS_tumor,WXS_normal,WXS_Status,RNA_seq_tumor,RNA_seq_Status)->SRP011540_three_paired



#SRP094781+SRP095809 WXS_normal~WXS_tumor~RNA_Seq
SRA %>%
  dplyr::filter(SRA_Study =="SRP094781")%>%
  dplyr::filter(Library_strategy == "RNA-Seq") %>%
  dplyr::select(Run,Biopsy_Time,patient_id) -> SRP094781
colnames(SRP094781)=c("RNA_seq_tumor","RNA_seq_Status","Patient_id")


paired_metadata %>%
  dplyr::filter(Project=="SRP095809")%>%
  merge(SRP094781,all=F)%>%
  dplyr::filter(WXS_Status==RNA_seq_Status)%>%
  dplyr::select(Cancer,Immune_checkpoint,Project,Patient_id,Response,Second_Response_standard,WXS_tumor,WXS_normal,WXS_Status,RNA_seq_tumor,RNA_seq_Status)->SRP094781_SRP095809_three_paired


#gastric cancer ERP107734 WXS_normal~WXS_tumor~RNA_Seq


#all
RNA_WXS_paired=rbind(SRP011540_three_paired,SRP094781_SRP095809_three_paired)
writexl::write_xlsx(RNA_WXS_paired,"/data/liucj/data/immune-checkpoint-blockade/RNA_WXS_paired.xlsx",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE,append=FALSE)
#caution about one timepoint two RNA-seq samples