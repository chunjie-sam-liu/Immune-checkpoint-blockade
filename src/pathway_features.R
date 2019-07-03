library(randomForest )


#features---------------------------------------------------------------------------
#mutant
DNA_damage=c("MSH2","MSH6","PMS2","POLE","BRCA2")
MAPK=c("KRAS","STK11","TP53")
JAK_STAT_mut=c("JAK1","JAK2","JAK3","SOCS1")
HLA_variability=c("B2M")
Chromatin_remodeling_mut=c("ARID1A","SMARCA4")

#protect
PI3K="PTEN"
PDL1="CD274"
JAK_STAT=c("IFNGR1","APLNR")
Urea_cycle_protect=c("SLC25A13","CPS1")
IFN_gamma=c("IFNG","STAT1","CXCL9","CXCL10","IDO1","HLA-DRA","CD3D","CIITA","CD3E","CCL5","GZMK","CD2",
            "CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","GZMB")
antigen_presentation=c("HLA-A","HLA-F","B2M","TAP1","TAP2")
T_cell_co_stimulation =c("ICAM1","CLECL1","LILRA1")#,"LILRA3"                   
cytokine=c("JAK2","STAT1")                     
cytotoxic_T=c("CD8A","CD8B","GZMA","GZMB","PRF1")


#risk
WNT_beta="DKK2"
IDO1 = "IDO1"
TGF_beta="TGFB1"
Chromatin_remodeling=c("PBRM1","EZH2")
Endo_retroviruses="KDM1A"
Urea_cycle_risk=c("ASL","OTC","SLC25A15","ASS1")

all_genes=c(PI3K,PDL1,JAK_STAT,Urea_cycle_protect,IFN_gamma,antigen_presentation,T_cell_co_stimulation,cytokine,cytotoxic_T,
            WNT_beta,IDO1,Chromatin_remodeling,Endo_retroviruses,Urea_cycle_risk,TGF_beta)
#load expr and metadata--------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_FPKM_expr.txt",
           header = T)->melanoma_PD1_FPKM
readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "SRA")%>%
  dplyr::filter(Cancer == "melanoma")%>%
  dplyr::filter(Anti_target =="anti-PD1")%>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(Run,SRA_Study,Response)%>%
  dplyr::filter(Response != "NE")->metadata
metadata$Response %>%
  gsub("PD","NR",.)%>% gsub("SD","NR",.)%>%
  gsub("PR","R",.)%>%gsub("CR","R",.)->metadata$Response


for (i in 1:ncol(melanoma_PD1_FPKM)) {
  which(is.na(melanoma_PD1_FPKM[,i]))->ID0
  melanoma_PD1_FPKM[ID0,i] = 0
}
sapply(melanoma_PD1_FPKM, function(x) log2(x+1)) %>%
  as.data.frame(stringsAsFactors=FALSE)->log2_melanoma_PD1_FPKM
rownames(log2_melanoma_PD1_FPKM)=rownames(melanoma_PD1_FPKM)

#1
metadata %>% dplyr::filter(SRA_Study=="SRP094781")->SRP094781_metadata
tibble::rownames_to_column(log2_melanoma_PD1_FPKM)%>%
  dplyr::select(c("rowname",SRP094781_metadata$Run))%>%
  tibble::column_to_rownames()->SRP094781_expr
scale(SRP094781_expr,center = TRUE,scale = TRUE)%>%as.data.frame(stringsAsFactors=FALSE)->scaled_SRP094781_expr

tibble::rownames_to_column(scaled_SRP094781_expr)%>%
  dplyr::filter(rowname %in% all_genes)%>%
  tibble::column_to_rownames()%>%
  t()%>%
  as.data.frame(stringsAsFactors=FALSE)->feature_expr

feature_score <- matrix(nrow = nrow(feature_expr),ncol = 14)
colnames(feature_score) = c("PI3K","PDL1","JAK_STAT","Urea_cycle","IFN_gamma","antigen_presentation","T_cell_co_stimulation",
                            "cytokine","cytotoxic_T","WNT_beta","IDO1","TGF_beta","Chromatin_remodeling","Endo_retroviruses")
rownames(feature_score)=rownames(feature_expr)
feature_score <- as.data.frame(feature_score,stringsAsFactors=FALSE)

feature_expr %>%
  dplyr::select(PI3K) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric()->feature_score$PI3K

feature_expr %>%
  dplyr::select(PDL1) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric()->feature_score$PDL1

feature_expr %>%
  dplyr::select(JAK_STAT) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric()->feature_score$JAK_STAT 
   
feature_expr %>%
  dplyr::select(IFN_gamma) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric()->feature_score$IFN_gamma  
  
feature_expr %>%
  dplyr::select(antigen_presentation) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric()->feature_score$antigen_presentation 

feature_expr %>%
  dplyr::select(T_cell_co_stimulation) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric()->feature_score$T_cell_co_stimulation 

feature_expr %>%
  dplyr::select(cytokine) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric()->feature_score$cytokine 

feature_expr %>%
  dplyr::select(cytotoxic_T) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric()->feature_score$cytotoxic_T 

feature_expr %>%
  dplyr::select(WNT_beta) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric() ->feature_score$WNT_beta 
 
feature_expr %>%
  dplyr::select(IDO1) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric() ->feature_score$IDO1

feature_expr %>%
  dplyr::select(TGF_beta) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric()->feature_score$TGF_beta

feature_expr %>%
  dplyr::select(Chromatin_remodeling) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric()->feature_score$Chromatin_remodeling

feature_expr %>%
  dplyr::select(Endo_retroviruses) %>%
  apply(., 1, mean) %>%
  as.character()%>%as.numeric()->feature_score$Endo_retroviruses

feature_expr %>%
  dplyr::select(Urea_cycle_protect,Urea_cycle_risk)%>%
  dplyr::mutate(Scores=(SLC25A13 + CPS1 - ASL - OTC - SLC25A15 - ASS1)/6)%>%
  dplyr::select(Scores)%>%
  as.matrix()%>%
  as.character()%>%as.numeric()-> feature_score$Urea_cycle


tibble::rownames_to_column(feature_score)%>%
  merge(SRP094781_metadata,by.x="rowname",by.y="Run")%>%
  dplyr::select(-SRA_Study)%>%
  tibble::column_to_rownames()->all_features_score

as.factor(all_features_score$Response)->all_features_score$Response

randomForest(x=all_features_score[,-15], y=all_features_score$Response, ntree = 1000,type="classification",
             importance=T )->features_rf



