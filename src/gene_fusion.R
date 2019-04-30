#PD1--------------------------------------------------------------------------------------

# path
path_data <- '/data/liull/immune-checkpoint-blockade/gene_fusion/SRP094781/'
# load gene fusion filtered.tsv
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(SRA_Study=="SRP094781") %>%
  dplyr::filter(Library_strategy=="RNA-Seq")%>%
  dplyr::filter(Anti_target=="anti-PD1")%>%#****
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata
dplyr::filter(metadata,Response %in% c("CR","PR","PRCR","R")) -> response#
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response#

#i=1

for (i in 1:nrow(response)) {
  list.files(path = path_data, pattern = response$Run[i], recursive = TRUE, full.names = TRUE)%>% 
    readr::read_tsv()  %>% 
    dplyr::select(GeneName1,GeneName2,fusiontype,GeneExpr1,GeneExpr2,GeneExpr_Fused,EricScore)-> filtered_pair
  apply(filtered_pair, 1, function(x) paste(x[1],x[2], sep = "_", collapse = NULL))->key
  dplyr::mutate(filtered_pair,key=key)->filtered_pair
  #print(nrow(filtered_pair))
  if(i==1){
    intersection=filtered_pair$key
  }else{
    intersection=intersect(filtered_pair$key,intersection) 
  }
  
}
#ERP107734:response:"C9orf172_PHPT1","STAG3L5P-PVRIG2P-PILRB_STAG3","ZNF547_TRAPPC2","RMND5A_ANAPC1","BPTF_KPNA2" 
#          nonresponse:"STAG3L5P-PVRIG2P-PILRB_STAG3"


#SRP070710:response:0
#          nonresponse:"C9orf172_PHPT1","STAG3L5P-PVRIG2P-PILRB_STAG3","ZNF547_TRAPPC2","PLEKHO2_ANKDD1A","BOLA2_SMG1"
#SRP150548:response:"ZNF547_TRAPPC2","SLC7A5_SMG1","DHRS1_THOP1","STAG3L5P-PVRIG2P-PILRB_STAG3","FBXO25_SEPT14","BPTF_KPNA2","EMC3_CIDEC","STYXL1_TMEM120A","SIDT2_TAGLN","PARG_TIMM23B","COX17_POPDC2"    
#          nonresponse:0
#SRP094781:response:KIAA0586_HNRNPC
#          nonresponse:0

#CTLA4----------------------------------------------------------------------------------------------------------------------------
# path
path_data <- '/data/liull/immune-checkpoint-blockade/gene_fusion/SRP011540/'
# load gene fusion filtered.tsv
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP") %>%
  dplyr::filter(SRA_Study=="SRP011540") %>%
  dplyr::filter(Library_strategy=="RNA-Seq")%>%
  dplyr::filter(Anti_target=="anti-CTLA4")%>%#****
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata
dplyr::filter(metadata,Response %in% c("CR","PR","X")) -> response#
dplyr::filter(metadata,Response %in% c("SD","PD")) -> non_response#

for (i in 1:nrow(non_response)) {
  list.files(path = path_data, pattern = non_response$Run[i], recursive = TRUE, full.names = TRUE)%>% 
    readr::read_tsv()  %>% 
    dplyr::select(GeneName1,GeneName2,fusiontype,GeneExpr1,GeneExpr2,GeneExpr_Fused,EricScore)-> filtered_pair
  apply(filtered_pair, 1, function(x) paste(x[1],x[2], sep = "_", collapse = NULL))->key
  dplyr::mutate(filtered_pair,key=key)->filtered_pair
  #print(nrow(filtered_pair))
  if(i==1){
    intersection=filtered_pair$key
  }else{
    intersection=intersect(filtered_pair$key,intersection) 
  }
  
}


#SRP011540:response:"STAG3L5P-PVRIG2P-PILRB_STAG3","HIC2_PI4KA","SLC7A5_SMG1" 
#          nonresponse:"STAG3L5P-PVRIG2P-PILRB_STAG3"
