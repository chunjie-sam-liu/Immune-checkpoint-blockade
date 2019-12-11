
options(stringsAsFactors=FALSE)
Projects=c("ERP107734","SRP070710","SRP094781","SRP150548","SRP011540")

#all metadata------------------------------------------------------------------
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP") %>%
  dplyr::filter(Cancer=="melanoma")%>%
  dplyr::filter(Anti_target=="anti-CTLA4")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(Run != "SRR3083584") %>%
  dplyr::select(SRA_Study,Run,Second_Response_standard) -> dbGAP
colnames(dbGAP)[3] <- "Response"
gsub("long-survival","R",dbGAP$Response) -> dbGAP$Response

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  #dplyr::filter(Cancer=="melanoma")%>%
  dplyr::filter(Anti_target=="anti-PD1")%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response)%>%
  dplyr::filter(Response != "NE")->SRA
gsub("PR","R",SRA$Response)%>%
  gsub("CR","R",.)%>%
  gsub("SD","NR",.)%>%
  gsub("PD","NR",.) -> SRA$Response
rbind(SRA,dbGAP)-> metadata


dplyr::filter(SRA_Study == "SRP094781")

for (i in 1:length(Projects)) {
  
  dplyr::filter(metadata,SRA_Study == Projects[i])->.metadata
  read.table(paste("/data/liull/immune-checkpoint-blockade/expression/",Projects[i],"_TIL.txt",sep = ""),
             header = T,as.is = TRUE) ->.TIL_file
  dplyr::filter(.metadata,Response == "R")->.R_ID
  dplyr::filter(.metadata,Response == "NR")->.NR_ID
  
  t(.TIL_file)%>%
    as.data.frame()%>%
    dplyr::select(.R_ID$Run,.NR_ID$Run)->.arranged_TIL
  
  P_value <- numeric()
  FC <- numeric()
  for (j in 1:nrow(.arranged_TIL)) {
    .arranged_TIL[j,1:nrow(.R_ID)]%>%as.numeric() -> R_TIL
    .arranged_TIL[j,(nrow(.R_ID)+1):nrow(.metadata)]%>%as.numeric() ->NR_TIL
    
    wilcox.test(R_TIL,NR_TIL,alternative = "two.sided")$p.value->.P_value
    P_value <- c(P_value,.P_value)
    sum(R_TIL)/sum(NR_TIL) ->.FC
    FC <- c(FC,.FC)
  }
  
  res <- data.frame(P_value = P_value,FC = FC)
  rownames(res) <- rownames(.arranged_TIL)
  write.table(res,paste("/data/liull/immune-checkpoint-blockade/TIL/TIL_R_VS_NR/",Projects[i],"_diff_res.txt",sep = ""),
              row.names = T,col.names = T,sep = "\t",quote = F)
  
}
  
  
