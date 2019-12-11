library(randomForest )
library()

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
Urea_cycle=c("SLC25A13","CPS1","ASL","OTC","SLC25A15","ASS1")#1,1,-1,-1,-1,-1
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


#genes
c("PTEN","CD274","IFNGR1","APLNR","SLC25A13","CPS1","ASL","OTC","SLC25A15","ASS1","IFNG","STAT1","CXCL9","CXCL10",
  "IDO1","HLA-DRA","CD3D","CIITA","CD3E","CCL5","GZMK","CD2","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP",
  "GZMB","HLA-A","HLA-F","B2M","TAP1","TAP2","ICAM1","CLECL1","LILRA1","JAK2","STAT1","CD8A","CD8B","GZMA","GZMB","PRF1",
  "DKK2","IDO1","TGFB1","PBRM1","EZH2","KDM1A")->signatures
#pairs
pairs=c("CD274","VSIR","CD40","PDCD1","CD28","CD80","TNFSF9",
        "CD86","CD276","TNFSF4","CD27","CD200","HAVCR2","TNFRSF14","CTLA4")
#VSIR = C10orf54

all_genes=c(PI3K,PDL1,JAK_STAT,Urea_cycle,IFN_gamma,antigen_presentation,T_cell_co_stimulation,cytokine,cytotoxic_T,
            WNT_beta,IDO1,TGF_beta,Chromatin_remodeling,Endo_retroviruses,pairs)


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
# sapply(melanoma_PD1_FPKM, function(x) log2(x+1)) %>%
#   as.data.frame(stringsAsFactors=FALSE)->log2_melanoma_PD1_FPKM
# rownames(log2_melanoma_PD1_FPKM)=rownames(melanoma_PD1_FPKM)

#train------------------------------------
metadata %>% dplyr::filter(SRA_Study=="SRP094781")->SRP094781_metadata

dplyr::filter(SRP094781_metadata,Response %in% c("CR","PR","R"))->responses
dplyr::filter(SRP094781_metadata,Response %in% c("PD","SD","NR"))->non_responses

tibble::rownames_to_column(melanoma_PD1_FPKM)%>%
  dplyr::select(c("rowname",responses$Run,non_responses$Run))%>%
  dplyr::filter(rowname %in% c(signatures,pairs))%>%
  tibble::column_to_rownames()->SRP094781_expr
scale(SRP094781_expr,center = TRUE,scale = TRUE)%>%as.data.frame(stringsAsFactors=FALSE)->scaled_SRP094781_expr

t(scaled_SRP094781_expr)%>%
  as.data.frame(stringsAsFactors=FALSE)->feature_expr

feature_score <- matrix(nrow = nrow(feature_expr),ncol = 14)
colnames(feature_score) = c("PI3K","PDL1","JAK_STAT","Urea_cycle","IFN_gamma","antigen_presentation","T_cell_co_stimulation",
                            "cytokine","cytotoxic_T","WNT_beta","IDO1","TGF_beta","Chromatin_remodeling","Endo_retroviruses")
rownames(feature_score)=rownames(feature_expr)
feature_score <- as.data.frame(feature_score,stringsAsFactors=FALSE)

# feature_expr %>%
#   dplyr::select(signatures)-> feature_score


fn_get_score <- function(expr,genes,index) {
  # prepare score dataframe.

  tibble::rownames_to_column(expr) %>%
    dplyr::select(c("rowname",genes)) %>%
    tibble::column_to_rownames()->scores_expr

  final_sum=rep(0,nrow(scores_expr))
  for(i in 1:length(genes)){

    final_sum = final_sum + scores_expr[,i]*index[i]

  }
  final_sum/length(genes)

}
# fn_pairs_score <- function(expr,gene1,gene2) {
#   # prepare pairs'score 
#   
#   tibble::rownames_to_column(expr) %>%
#     dplyr::select(c("rowname",gene1,gene2)) %>%
#     tibble::column_to_rownames()->scores_expr
#   
#   res=rep(0,nrow(scores_expr))
#   for(i in 1:nrow(scores_expr)){
#     
#     if(scores_expr[i,1] > scores_expr[i,2]){
#       res[i]=1
#     }
#     
#   }
#   res
# }

fn_pairs_score <- function(expr,gene1,gene2) {
  # prepare pairs'score

  tibble::rownames_to_column(expr) %>%
    dplyr::select(c("rowname",gene1,gene2)) %>%
    tibble::column_to_rownames()->scores_expr

  res=rep(0,nrow(scores_expr))
  for(i in 1:nrow(scores_expr)){

    res[i]=scores_expr[i,1]/scores_expr[i,2]

  }
  res
}

feature_score$PI3K=fn_get_score(feature_expr,PI3K,1)
feature_score$PDL1=fn_get_score(feature_expr,PDL1,1)
feature_score$JAK_STAT=fn_get_score(feature_expr,JAK_STAT,c(1,1))
feature_score$Urea_cycle=fn_get_score(feature_expr,Urea_cycle,c(1,1,-1,-1,-1,-1))
feature_score$IFN_gamma=fn_get_score(feature_expr,IFN_gamma,rep(1,20))
feature_score$antigen_presentation=fn_get_score(feature_expr,antigen_presentation,rep(1,5))
feature_score$T_cell_co_stimulation=fn_get_score(feature_expr,T_cell_co_stimulation,rep(1,3))
feature_score$cytokine=fn_get_score(feature_expr,cytokine,c(1,1))
feature_score$cytotoxic_T=fn_get_score(feature_expr,cytotoxic_T,rep(1,5))

feature_score$WNT_beta=fn_get_score(feature_expr,WNT_beta,-1)
feature_score$IDO1=fn_get_score(feature_expr,IDO1,-1)
feature_score$TGF_beta=fn_get_score(feature_expr,TGF_beta,-1)
feature_score$Chromatin_remodeling=fn_get_score(feature_expr,Chromatin_remodeling,rep(-1,2))
feature_score$Endo_retroviruses=fn_get_score(feature_expr,Endo_retroviruses,-1)

tibble::rownames_to_column(feature_score) %>%
  dplyr::mutate(CD274_VSIR=fn_pairs_score(feature_expr,"CD274","VSIR"))%>%
  dplyr::mutate(CD40_CD274=fn_pairs_score(feature_expr,"CD40","CD274"))%>%
  dplyr::mutate(CD40_PDCD1=fn_pairs_score(feature_expr,"CD40","PDCD1"))%>%
  dplyr::mutate(CD40_CD28=fn_pairs_score(feature_expr,"CD40","CD28"))%>%
  dplyr::mutate(CD40_CD80=fn_pairs_score(feature_expr,"CD40","CD80"))%>%
  dplyr::mutate(CD80_TNFSF9=fn_pairs_score(feature_expr,"CD80","TNFSF9"))%>%
  dplyr::mutate(CD28_CD86=fn_pairs_score(feature_expr,"CD28","CD86"))%>%
  dplyr::mutate(CD28_CD276=fn_pairs_score(feature_expr,"CD28","CD276"))%>%
  dplyr::mutate(PDCD1_TNFSF4=fn_pairs_score(feature_expr,"PDCD1","TNFSF4"))%>%
  dplyr::mutate(CD27_PDCD1=fn_pairs_score(feature_expr,"CD27","PDCD1"))%>%
  dplyr::mutate(CD86_TNFSF4=fn_pairs_score(feature_expr,"CD86","TNFSF4"))%>%
  dplyr::mutate(CD86_CD200=fn_pairs_score(feature_expr,"CD86","CD200"))%>%
  dplyr::mutate(CD86_HAVCR2=fn_pairs_score(feature_expr,"CD86","HAVCR2"))%>%
  dplyr::mutate(TNFRSF14_CD86=fn_pairs_score(feature_expr,"TNFRSF14","CD86"))%>%
  dplyr::mutate(CTLA4_TNFSF4=fn_pairs_score(feature_expr,"CTLA4","TNFSF4"))->feature_score



# read.table("/data/liull/immune-checkpoint-blockade/expression/SRP094781_TIL.txt",header = T)->SRP094781_TIL
# tibble::rownames_to_column(SRP094781_TIL)%>%
#   merge(feature_score,by="rowname")->feature_score


merge(SRP094781_metadata,feature_score,by.x="Run",by.y="rowname")%>%
  dplyr::select(-SRA_Study)->all_features_score
all_features_score$Run -> rownames(all_features_score)
all_features_score %>% dplyr::select(-Run) ->all_features_score

as.factor(all_features_score$Response)->all_features_score$Response


#randomForest
randomForest(x=all_features_score[,-1], y=all_features_score$Response, ntree = 500,type="classification",
             importance=T )->features_rf

#svm
svm(Response ~ ., data = all_features_score,type = 'C',kernel = 'radial' )->features_svm
pred_svm <- predict(features_svm,all_features_score,type="prob")
dim(all_features_score[all_features_score$Response!=pred_svm,])[1]/dim(all_features_score)[1]#error rate


# train_model_rf=train(x=all_features_score[,-30],y=all_features_score$Response,
#                       trControl=trainControl(method = "repeatedcv",repeats = 5,classProbs = TRUE),
#                       method="rf",family=binomial, metric = "ROC")
# pred_train=predict(train_model_rf,all_features_score[,-30],type="prob")
# roc_train=roc(predictor=pred_train[,1],response=all_features_score$Response,levels=rev(levels(all_features_score$Response)),plot=T,ci=T,xlab="1 - Specificity")


#test---------------------------------------------------------------------------------
metadata %>% dplyr::filter(SRA_Study=="SRP070710")->SRP070710_metadata
tibble::rownames_to_column(melanoma_PD1_FPKM)%>%
  dplyr::select(c("rowname",SRP070710_metadata$Run))%>%
  tibble::column_to_rownames()->SRP070710_expr
scale(SRP070710_expr,center = TRUE,scale = TRUE)%>%as.data.frame(stringsAsFactors=FALSE)->scaled_SRP070710_expr

tibble::rownames_to_column(scaled_SRP070710_expr)%>%
  dplyr::filter(rowname %in% all_genes)%>%
  tibble::column_to_rownames()%>%
  t()%>%
  as.data.frame(stringsAsFactors=FALSE)->test_expr

test_score <- matrix(nrow = nrow(test_expr),ncol = 14)
colnames(test_score) = c("PI3K","PDL1","JAK_STAT","Urea_cycle","IFN_gamma","antigen_presentation","T_cell_co_stimulation",
                            "cytokine","cytotoxic_T","WNT_beta","IDO1","TGF_beta","Chromatin_remodeling","Endo_retroviruses")
rownames(test_score)=rownames(test_expr)
test_score <- as.data.frame(test_score,stringsAsFactors=FALSE)

test_score$PI3K=fn_get_score(test_expr,PI3K,1)
test_score$PDL1=fn_get_score(test_expr,PDL1,1)
test_score$JAK_STAT=fn_get_score(test_expr,JAK_STAT,c(1,1))
test_score$Urea_cycle=fn_get_score(test_expr,Urea_cycle,c(1,1,-1,-1,-1,-1))
test_score$IFN_gamma=fn_get_score(test_expr,IFN_gamma,rep(1,20))
test_score$antigen_presentation=fn_get_score(test_expr,antigen_presentation,rep(1,5))
test_score$T_cell_co_stimulation=fn_get_score(test_expr,T_cell_co_stimulation,rep(1,3))
test_score$cytokine=fn_get_score(test_expr,cytokine,c(1,1))
test_score$cytotoxic_T=fn_get_score(test_expr,cytotoxic_T,rep(1,5))

test_score$WNT_beta=fn_get_score(test_expr,WNT_beta,-1)
test_score$IDO1=fn_get_score(test_expr,IDO1,-1)
test_score$TGF_beta=fn_get_score(test_expr,TGF_beta,-1)
test_score$Chromatin_remodeling=fn_get_score(test_expr,Chromatin_remodeling,rep(-1,2))
test_score$Endo_retroviruses=fn_get_score(test_expr,Endo_retroviruses,-1)


tibble::rownames_to_column(test_score) %>%
  dplyr::mutate(pair1=fn_pairs_score(test_expr,"CD274","VSIR"))%>%
  dplyr::mutate(pair2=fn_pairs_score(test_expr,"CD40","CD274"))%>%
  dplyr::mutate(pair3=fn_pairs_score(test_expr,"CD40","PDCD1"))%>%
  dplyr::mutate(pair4=fn_pairs_score(test_expr,"CD40","CD28"))%>%
  dplyr::mutate(pair5=fn_pairs_score(test_expr,"CD40","CD80"))%>%
  dplyr::mutate(pair6=fn_pairs_score(test_expr,"CD80","TNFSF9"))%>%
  dplyr::mutate(pair7=fn_pairs_score(test_expr,"CD28","CD86"))%>%
  dplyr::mutate(pair8=fn_pairs_score(test_expr,"CD28","CD276"))%>%
  dplyr::mutate(pair9=fn_pairs_score(test_expr,"PDCD1","TNFSF4"))%>%
  dplyr::mutate(pair10=fn_pairs_score(test_expr,"CD27","PDCD1"))%>%
  dplyr::mutate(pair11=fn_pairs_score(test_expr,"CD86","TNFSF4"))%>%
  dplyr::mutate(pair12=fn_pairs_score(test_expr,"CD86","CD200"))%>%
  dplyr::mutate(pair13=fn_pairs_score(test_expr,"CD86","HAVCR2"))%>%
  dplyr::mutate(pair14=fn_pairs_score(test_expr,"TNFRSF14","CD86"))%>%
  dplyr::mutate(pair15=fn_pairs_score(test_expr,"CTLA4","TNFSF4"))->test_score

read.table("/data/liull/immune-checkpoint-blockade/expression/SRP070710_TIL.txt",header = T)->SRP070710_TIL
tibble::rownames_to_column(SRP070710_TIL)%>%
  merge(test_score,by="rowname")->test_score


merge(SRP070710_metadata,test_score,by.x="Run",by.y="rowname")%>%
  dplyr::select(-SRA_Study)->all_test_score
all_test_score$Run -> rownames(all_test_score)
all_test_score %>% dplyr::select(-Run) ->all_test_score
as.factor(all_test_score$Response)->all_test_score$Response

#randomForst
pred_test=predict(features_rf,all_test_score[,-1])
ran_roc <- roc(all_test_score$Response,as.numeric(pred_test))
plot(ran_roc, print.auc=TRUE, auc.polygon=TRUE,max.auc.polygon=TRUE,auc.polygon.col="skyblue")

#svm
pred_test=predict(features_svm,all_test_score)
ran_roc <- roc(all_test_score$Response,as.numeric(pred_test))
plot(ran_roc, print.auc=TRUE, auc.polygon=TRUE,max.auc.polygon=TRUE,auc.polygon.col="skyblue")

# pred_test=predict(features_svm,all_test_score,type="prob")
# dim(all_test_score[all_test_score$Response!=pred_test,])[1]/dim(all_test_score)[1]#error rate





#melanoam_PD1 3 project together--------------------------------------- 
library(GSVA)

#features

#protect
PI3K="PTEN"
PDL1="CD274"
JAK_STAT=c("IFNGR1","APLNR")
Urea_cycle=c("SLC25A13","CPS1","ASL","OTC","SLC25A15","ASS1")#1,1,-1,-1,-1,-1
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
#pairs
pairs=c()


all_genes=c("PTEN","CD274","IFNGR1","APLNR","SLC25A13","CPS1","ASL","OTC","SLC25A15","ASS1","IFNG","STAT1","CXCL9","CXCL10","IDO1","HLA-DRA",
            "CD3D","CIITA","CD3E","CCL5","GZMK","CD2","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","GZMB","HLA-A","HLA-F","B2M","TAP1",
            "TAP2","ICAM1","CLECL1","LILRA1","JAK2","STAT1","CD8A","CD8B","GZMA","GZMB","PRF1","DKK2","IDO1","TGFB1","PBRM1","EZH2","KDM1A",
            "CD274","VSIR","CD40","PDCD1","CD28","CD80","TNFSF9","CD86","CD276","TNFSF4", "CD27","CD200","HAVCR2","TNFRSF14","CTLA4")

gene_sets=list()
gene_sets$PI3K="PTEN"
gene_sets$PDL1="CD274"
gene_sets$JAK_STAT=c("IFNGR1","APLNR")
gene_sets$Urea_cycle=c("SLC25A13","CPS1","ASL","OTC","SLC25A15","ASS1")
gene_sets$IFN_gamma=c("IFNG","STAT1","CXCL9","CXCL10","IDO1","HLA-DRA","CD3D","CIITA","CD3E","CCL5","GZMK","CD2",
               "CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","GZMB")
gene_sets$antigen_presentation=c("HLA-A","HLA-F","B2M","TAP1","TAP2")
gene_sets$T_cell_co_stimulation=c("ICAM1","CLECL1","LILRA1")
gene_sets$cytokine=c("JAK2","STAT1")
gene_sets$cytotoxic_T=c("CD8A","CD8B","GZMA","GZMB","PRF1")
gene_sets$WNT_beta="DKK2"
gene_sets$IDO1="IDO1"
gene_sets$TGF_beta="TGFB1"
gene_sets$Chromatin_remodeling=c("PBRM1","EZH2")
gene_sets$Endo_retroviruses="KDM1A"

gene_sets$CD274_VSIR=c("CD274","VSIR")
gene_sets$CD40_CD274=c("CD40","CD274")
gene_sets$CD40_PDCD1=c("CD40","PDCD1")
gene_sets$CD40_CD28=c("CD40","CD28")
gene_sets$CD40_CD80=c("CD40","CD80")
gene_sets$CD80_TNFSF9=c("CD80","TNFSF9")
gene_sets$CD28_CD86=c("CD28","CD86")
gene_sets$CD28_CD276=c("CD28","CD276")
gene_sets$PDCD1_TNFSF4=c("PDCD1","TNFSF4")
gene_sets$CD27_PDCD1=c("CD27","PDCD1")
gene_sets$CD86_TNFSF4=c("CD86","TNFSF4")
gene_sets$CD86_CD200=c("CD86","CD200")
gene_sets$CD86_HAVCR2=c("CD86","HAVCR2")
gene_sets$TNFRSF14_CD86=c("TNFRSF14","CD86")
gene_sets$CTLA4_TNFSF4=c("CTLA4","TNFSF4")


#load expr to ssgsva and test in response and nonresponse--------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_log2CPM_expr.txt",
           header = T)->melanoma_PD1_log2CPM

ssgava_score <- gsva(as.matrix(melanoma_PD1_log2CPM), gene_sets, min.sz=1, max.sz=999999, method="ssgsea",kcdf="Gaussian")


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


#boxplot--------------------------------------------------------------------------------------
t(ssgava_score)%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  tibble::rownames_to_column()%>%
  merge(metadata,.,by.x="Run",by.y="rowname")%>%
  dplyr::select(-SRA_Study,-Run)->plot_score



for (j in 1:5) {
  
  file_names=paste("/data/liull/immune-checkpoint-blockade/pathway_features/features_",j,".pdf",sep = "")
  
  pdf(file_names)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3,2)))
  
  vplayout <- function(x,y){
    viewport(layout.pos.row = x, layout.pos.col = y)
  }
  
  for (i in (6*j-5):(6*j)) {   #while j=5,x=25:29
    
    if(i%%2==0){
      n=2
    }else{
      n=1
    }  
    
    if(i%%6==1 || i%%6==2){
      m=1
    }else if(i%%6==3 || i%%6==4){
      m=2
    }else{
      m=3
    }
    
    dplyr::select(plot_score,Response,colnames(plot_score)[i+1])->single_feature_expr
    
    dplyr::filter(single_feature_expr,Response=="R") %>% dplyr::select(-Response)->R_expr
    dplyr::filter(single_feature_expr,Response=="NR") %>% dplyr::select(-Response)->NR_expr
    
    wilcox.test(R_expr[,1],NR_expr[,1])$p.value -> anno
    
    print(
      ggplot(single_feature_expr)+
        geom_boxplot(aes(x=Response,y=as.numeric(single_feature_expr[,2]),fill=Response))+
        labs(y=colnames(plot_score)[i+1])+
        annotate(geom="text", x=1.5, y=(max(as.numeric(single_feature_expr[,2])))/2,size=3, label=paste("P.value",signif(anno,2),sep = ":"))+
        scale_fill_manual(values=c("#EE7600","#FFD700"))+
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border=element_rect(fill = NA),
              axis.text.x = element_text(size = 5)),
      vp = vplayout(m,n)
    )
    
  }
  dev.off()
}

#test our data's problem with SRP070710 in GEO------------------
readxl::read_xlsx("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/GSE78220_PatientFPKM.xlsx")->SRP070710_GEO

readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "SRA")%>%
  dplyr::filter(Second_accession == "GSE78220")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(Run,Response,patient_id)->metadata
gsub(".baseline","",colnames(SRP070710_GEO))%>%gsub(".OnTx","",.)%>%
  gsub("A","",.)%>%gsub("B","",.)->colnames(SRP070710_GEO)

for (i in 1:nrow(metadata)) {
  
  which(colnames(SRP070710_GEO)[] == metadata$patient_id[i])->ID
  colnames(SRP070710_GEO)[ID] <- metadata$Run[i]
}
colnames(SRP070710_GEO)[22]="SRR3184299"

dplyr::select(SRP070710_GEO,"Gene",metadata$Run)%>%
  as.data.frame(stringsAsFactors=FALSE)->SRP070710_GEO
rownames(SRP070710_GEO) <- SRP070710_GEO$Gene
SRP070710_GEO[,-1]->SRP070710_GEO


#zscore by fn_------------------------------------------
scale(SRP070710_GEO,center = TRUE,scale = TRUE)%>%
  t()%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  dplyr::select(all_genes)->scaled_SRP070710_GEO



GEO_feature_score <- matrix(nrow = nrow(scaled_SRP070710_GEO),ncol = 14)
colnames(GEO_feature_score) = c("PI3K","PDL1","JAK_STAT","Urea_cycle","IFN_gamma","antigen_presentation","T_cell_co_stimulation",
                                "cytokine","cytotoxic_T","WNT_beta","IDO1","TGF_beta","Chromatin_remodeling","Endo_retroviruses")
rownames(GEO_feature_score)=rownames(scaled_SRP070710_GEO)
GEO_feature_score <- as.data.frame(GEO_feature_score,stringsAsFactors=FALSE)


GEO_feature_score$PI3K=fn_get_score(scaled_SRP070710_GEO,PI3K,1)
GEO_feature_score$PDL1=fn_get_score(scaled_SRP070710_GEO,PDL1,1)
GEO_feature_score$JAK_STAT=fn_get_score(scaled_SRP070710_GEO,JAK_STAT,c(1,1))
GEO_feature_score$Urea_cycle=fn_get_score(scaled_SRP070710_GEO,Urea_cycle,c(1,1,-1,-1,-1,-1))
GEO_feature_score$IFN_gamma=fn_get_score(scaled_SRP070710_GEO,IFN_gamma,rep(1,20))
GEO_feature_score$antigen_presentation=fn_get_score(scaled_SRP070710_GEO,antigen_presentation,rep(1,5))
GEO_feature_score$T_cell_co_stimulation=fn_get_score(scaled_SRP070710_GEO,T_cell_co_stimulation,rep(1,3))
GEO_feature_score$cytokine=fn_get_score(scaled_SRP070710_GEO,cytokine,c(1,1))
GEO_feature_score$cytotoxic_T=fn_get_score(scaled_SRP070710_GEO,cytotoxic_T,rep(1,5))

GEO_feature_score$WNT_beta=fn_get_score(scaled_SRP070710_GEO,WNT_beta,-1)
GEO_feature_score$IDO1=fn_get_score(scaled_SRP070710_GEO,IDO1,-1)
GEO_feature_score$TGF_beta=fn_get_score(scaled_SRP070710_GEO,TGF_beta,-1)
GEO_feature_score$Chromatin_remodeling=fn_get_score(scaled_SRP070710_GEO,Chromatin_remodeling,rep(-1,2))
GEO_feature_score$Endo_retroviruses=fn_get_score(scaled_SRP070710_GEO,Endo_retroviruses,-1)

tibble::rownames_to_column(GEO_feature_score) %>%
  dplyr::mutate(CD274_VSIR=fn_pairs_score(scaled_SRP070710_GEO,"CD274","C10orf54"))%>%
  dplyr::mutate(CD40_CD274=fn_pairs_score(scaled_SRP070710_GEO,"CD40","CD274"))%>%
  dplyr::mutate(CD40_PDCD1=fn_pairs_score(scaled_SRP070710_GEO,"CD40","PDCD1"))%>%
  dplyr::mutate(CD40_CD28=fn_pairs_score(scaled_SRP070710_GEO,"CD40","CD28"))%>%
  dplyr::mutate(CD40_CD80=fn_pairs_score(scaled_SRP070710_GEO,"CD40","CD80"))%>%
  dplyr::mutate(CD80_TNFSF9=fn_pairs_score(scaled_SRP070710_GEO,"CD80","TNFSF9"))%>%
  dplyr::mutate(CD28_CD86=fn_pairs_score(scaled_SRP070710_GEO,"CD28","CD86"))%>%
  dplyr::mutate(CD28_CD276=fn_pairs_score(scaled_SRP070710_GEO,"CD28","CD276"))%>%
  dplyr::mutate(PDCD1_TNFSF4=fn_pairs_score(scaled_SRP070710_GEO,"PDCD1","TNFSF4"))%>%
  dplyr::mutate(CD27_PDCD1=fn_pairs_score(scaled_SRP070710_GEO,"CD27","PDCD1"))%>%
  dplyr::mutate(CD86_TNFSF4=fn_pairs_score(scaled_SRP070710_GEO,"CD86","TNFSF4"))%>%
  dplyr::mutate(CD86_CD200=fn_pairs_score(scaled_SRP070710_GEO,"CD86","CD200"))%>%
  dplyr::mutate(CD86_HAVCR2=fn_pairs_score(scaled_SRP070710_GEO,"CD86","HAVCR2"))%>%
  dplyr::mutate(TNFRSF14_CD86=fn_pairs_score(scaled_SRP070710_GEO,"TNFRSF14","CD86"))%>%
  dplyr::mutate(CTLA4_TNFSF4=fn_pairs_score(scaled_SRP070710_GEO,"CTLA4","TNFSF4"))->GEO_feature_score

metadata$Response %>%
  gsub("PD","NR",.)%>% gsub("SD","NR",.)%>%
  gsub("PR","R",.)%>%gsub("CR","R",.)->metadata$Response
merge(GEO_feature_score,metadata,by.x="rowname",by.y="Run")%>%
  dplyr::select(-patient_id)%>%
  tibble::column_to_rownames()->GEO_score

for (i in 1:6) {
  
  dplyr::select(GEO_score,Response,colnames(GEO_score)[i])->single_feature_expr
  
  dplyr::filter(single_feature_expr,Response == "R") %>% dplyr::select(-Response)->R_expr
  dplyr::filter(single_feature_expr,Response == "NR") %>% dplyr::select(-Response)->NR_expr
  
  wilcox.test(R_expr[,1],NR_expr[,1])$p.value -> anno
  print(anno)
  ggplot(single_feature_expr)+
    geom_boxplot(aes(x=Response,y=as.numeric(single_feature_expr[,2]),fill=Response))+
    labs(y=colnames(GEO_score)[i+1])+
    annotate(geom="text", x=1.5, y=(max(as.numeric(single_feature_expr[,2])))/2,size=3, label=paste("P.value",signif(anno,2),sep = ":"))+
    scale_fill_manual(values=c("#EE7600","#FFD700"))+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border=element_rect(fill = NA),
          axis.text.x = element_text(size = 5))
  
}

##zscore by GSVA------------------------------------------
all_genes=c("PTEN","CD274","IFNGR1","APLNR","SLC25A13","CPS1","ASL","OTC","SLC25A15","ASS1","IFNG","STAT1","CXCL9","CXCL10","IDO1","HLA-DRA",
            "CD3D","CIITA","CD3E","CCL5","GZMK","CD2","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","GZMB","HLA-A","HLA-F","B2M","TAP1",
            "TAP2","ICAM1","CLECL1","LILRA1","JAK2","STAT1","CD8A","CD8B","GZMA","GZMB","PRF1","DKK2","IDO1","TGFB1","PBRM1","EZH2","KDM1A",
            "CD274","C10orf54","CD40","PDCD1","CD28","CD80","TNFSF9","CD86","CD276","TNFSF4", "CD27","CD200","HAVCR2","TNFRSF14","CTLA4")

gene_sets=list()
gene_sets$PI3K="PTEN"
gene_sets$PDL1="CD274"
gene_sets$JAK_STAT=c("IFNGR1","APLNR")
gene_sets$Urea_cycle=c("SLC25A13","CPS1","ASL","OTC","SLC25A15","ASS1")
gene_sets$IFN_gamma=c("IFNG","STAT1","CXCL9","CXCL10","IDO1","HLA-DRA","CD3D","CIITA","CD3E","CCL5","GZMK","CD2",
                      "CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","GZMB")
gene_sets$antigen_presentation=c("HLA-A","HLA-F","B2M","TAP1","TAP2")
gene_sets$T_cell_co_stimulation=c("ICAM1","CLECL1","LILRA1")
gene_sets$cytokine=c("JAK2","STAT1")
gene_sets$cytotoxic_T=c("CD8A","CD8B","GZMA","GZMB","PRF1")
gene_sets$WNT_beta="DKK2"
gene_sets$IDO1="IDO1"
gene_sets$TGF_beta="TGFB1"
gene_sets$Chromatin_remodeling=c("PBRM1","EZH2")
gene_sets$Endo_retroviruses="KDM1A"

gene_sets$CD274_VSIR=c("CD274","C10orf54")
gene_sets$CD40_CD274=c("CD40","CD274")
gene_sets$CD40_PDCD1=c("CD40","PDCD1")
gene_sets$CD40_CD28=c("CD40","CD28")
gene_sets$CD40_CD80=c("CD40","CD80")
gene_sets$CD80_TNFSF9=c("CD80","TNFSF9")
gene_sets$CD28_CD86=c("CD28","CD86")
gene_sets$CD28_CD276=c("CD28","CD276")
gene_sets$PDCD1_TNFSF4=c("PDCD1","TNFSF4")
gene_sets$CD27_PDCD1=c("CD27","PDCD1")
gene_sets$CD86_TNFSF4=c("CD86","TNFSF4")
gene_sets$CD86_CD200=c("CD86","CD200")
gene_sets$CD86_HAVCR2=c("CD86","HAVCR2")
gene_sets$TNFRSF14_CD86=c("TNFRSF14","CD86")
gene_sets$CTLA4_TNFSF4=c("CTLA4","TNFSF4")

readxl::read_xlsx("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/GSE78220_PatientFPKM.xlsx")->SRP070710_GEO

readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "SRA")%>%
  dplyr::filter(Second_accession == "GSE78220")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(Run,Response,patient_id)->metadata
gsub(".baseline","",colnames(SRP070710_GEO))%>%gsub(".OnTx","",.)%>%
  gsub("A","",.)%>%gsub("B","",.)->colnames(SRP070710_GEO)

for (i in 1:nrow(metadata)) {
  
  which(colnames(SRP070710_GEO)[] == metadata$patient_id[i])->ID
  colnames(SRP070710_GEO)[ID] <- metadata$Run[i]
}
colnames(SRP070710_GEO)[22]="SRR3184299"

dplyr::select(SRP070710_GEO,"Gene",metadata$Run)%>%
  as.data.frame(stringsAsFactors=FALSE)->SRP070710_GEO
rownames(SRP070710_GEO) <- SRP070710_GEO$Gene
SRP070710_GEO[,-1]->SRP070710_GEO

gava_zscore <- gsva(as.matrix(SRP070710_GEO), gene_sets, min.sz=1, max.sz=999999, method="zscore",kcdf="Gaussian")

metadata$Response %>%
  gsub("PD","NR",.)%>% gsub("SD","NR",.)%>%
  gsub("PR","R",.)%>%gsub("CR","R",.)->metadata$Response

t(gava_zscore)%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  tibble::rownames_to_column()%>%
  merge(metadata,.,by.x="Run",by.y="rowname")%>%
  dplyr::select(-patient_id,-Run)->plot_score



for (j in 1:5) {
  
  file_names=paste("/data/liull/immune-checkpoint-blockade/pathway_features/features_",j,".pdf",sep = "")
  
  pdf(file_names)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3,2)))
  
  vplayout <- function(x,y){
    viewport(layout.pos.row = x, layout.pos.col = y)
  }
  
  for (i in (6*j-5):(6*j)) {   #while j=5,x=25:29
    
    if(i%%2==0){
      n=2
    }else{
      n=1
    }  
    
    if(i%%6==1 || i%%6==2){
      m=1
    }else if(i%%6==3 || i%%6==4){
      m=2
    }else{
      m=3
    }
    
    dplyr::select(plot_score,Response,colnames(plot_score)[i+1])->single_feature_expr
    
    dplyr::filter(single_feature_expr,Response=="R") %>% dplyr::select(-Response)->R_expr
    dplyr::filter(single_feature_expr,Response=="NR") %>% dplyr::select(-Response)->NR_expr
    
    wilcox.test(R_expr[,1],NR_expr[,1])$p.value -> anno
    
    print(
      ggplot(single_feature_expr)+
        geom_boxplot(aes(x=Response,y=as.numeric(single_feature_expr[,2]),fill=Response))+
        labs(y=colnames(plot_score)[i+1])+
        annotate(geom="text", x=1.5, y=(max(as.numeric(single_feature_expr[,2])))/2,size=3, label=paste("P.value",signif(anno,2),sep = ":"))+
        scale_fill_manual(values=c("#EE7600","#FFD700"))+
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border=element_rect(fill = NA),
              axis.text.x = element_text(size = 5)),
      vp = vplayout(m,n)
    )
    
  }
  dev.off()
}

