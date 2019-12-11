fn_get_score <- function(expr,genes,index) {
  # prepare score dataframe.
  intersect(genes,colnames(expr)) ->genes
  
  if(length(genes) >0){
    tibble::rownames_to_column(expr) %>%
      dplyr::select(c("rowname",genes)) %>%
      tibble::column_to_rownames()->scores_expr
    
    final_sum=rep(0,nrow(scores_expr))
    for(i in 1:length(genes)){
      
      final_sum = final_sum + scores_expr[,i]*index[i]
    }
    res <- final_sum/length(genes)
  }else{
    res <- rep(0,nrow(expr))
  }
  res
}


fn_pairs_score <- function(expr,gene1,gene2) {
  # prepare pairs'score
  
  if(length(intersect(gene1,colnames(expr))) == 1 && length(intersect(gene2,colnames(expr))) == 1){
    
    tibble::rownames_to_column(expr) %>%
      dplyr::select(c("rowname",gene1,gene2)) %>%
      tibble::column_to_rownames()->scores_expr
    
    res=rep(0,nrow(scores_expr))
    for(i in 1:nrow(scores_expr)){
      
      res[i] <- scores_expr[i,1]/scores_expr[i,2]
      
    }
    
  }else if(length(intersect(gene1,colnames(expr))) == 1 && length(intersect(gene2,colnames(expr))) != 1){
    dplyr::select(expr,gene1) %>%
      as.matrix()%>%as.numeric()->res
  }else if(length(intersect(gene1,colnames(expr))) != 1 && length(intersect(gene2,colnames(expr))) == 1){
    dplyr::select(expr,gene2) %>%
      as.matrix()%>%as.numeric()->res
    1/res ->res
  }else{
    res <- rep(0,nrow(expr))
  }
  
  res
}

fn_arrange_score_matrix <- function(scaled_expr){
  
  Project_score <- matrix(nrow = nrow(scaled_expr),ncol = 14)%>%as.data.frame()
  colnames(Project_score) = c("PI3K","PDL1","JAK_STAT","Urea_cycle","IFN_gamma","antigen_presentation","T_cell_co_stimulation",
                              "cytokine","cytotoxic_T","WNT_beta","IDO1","TGF_beta","Chromatin_remodeling","Endo_retroviruses")
  rownames(Project_score)=rownames(scaled_expr)
  
  
  Project_score$PI3K=fn_get_score(scaled_expr,PI3K,1)
  Project_score$PDL1=fn_get_score(scaled_expr,PDL1,1)
  Project_score$JAK_STAT=fn_get_score(scaled_expr,JAK_STAT,c(1,1))
  Project_score$Urea_cycle=fn_get_score(scaled_expr,Urea_cycle,c(1,1,-1,-1,-1,-1))
  Project_score$IFN_gamma=fn_get_score(scaled_expr,IFN_gamma,rep(1,20))
  Project_score$antigen_presentation=fn_get_score(scaled_expr,antigen_presentation,rep(1,5))
  Project_score$T_cell_co_stimulation=fn_get_score(scaled_expr,T_cell_co_stimulation,rep(1,3))
  Project_score$cytokine=fn_get_score(scaled_expr,cytokine,c(1,1))
  Project_score$cytotoxic_T=fn_get_score(scaled_expr,cytotoxic_T,rep(1,5))
  
  Project_score$WNT_beta=fn_get_score(scaled_expr,WNT_beta,-1)
  Project_score$IDO1=fn_get_score(scaled_expr,IDO1,-1)
  Project_score$TGF_beta=fn_get_score(scaled_expr,TGF_beta,-1)
  Project_score$Chromatin_remodeling=fn_get_score(scaled_expr,Chromatin_remodeling,rep(-1,2))
  Project_score$Endo_retroviruses=fn_get_score(scaled_expr,Endo_retroviruses,-1)
  
  tibble::rownames_to_column(Project_score) %>%
    dplyr::mutate(CD274_VSIR=fn_pairs_score(scaled_expr,"CD274","VSIR"))%>%
    dplyr::mutate(CD40_CD274=fn_pairs_score(scaled_expr,"CD40","CD274"))%>%
    dplyr::mutate(CD40_PDCD1=fn_pairs_score(scaled_expr,"CD40","PDCD1"))%>%
    dplyr::mutate(CD40_CD28=fn_pairs_score(scaled_expr,"CD40","CD28"))%>%
    dplyr::mutate(CD40_CD80=fn_pairs_score(scaled_expr,"CD40","CD80"))%>%
    dplyr::mutate(CD80_TNFSF9=fn_pairs_score(scaled_expr,"CD80","TNFSF9"))%>%
    dplyr::mutate(CD28_CD86=fn_pairs_score(scaled_expr,"CD28","CD86"))%>%
    dplyr::mutate(CD28_CD276=fn_pairs_score(scaled_expr,"CD28","CD276"))%>%
    dplyr::mutate(PDCD1_TNFSF4=fn_pairs_score(scaled_expr,"PDCD1","TNFSF4"))%>%
    dplyr::mutate(CD27_PDCD1=fn_pairs_score(scaled_expr,"CD27","PDCD1"))%>%
    dplyr::mutate(CD86_TNFSF4=fn_pairs_score(scaled_expr,"CD86","TNFSF4"))%>%
    dplyr::mutate(CD86_CD200=fn_pairs_score(scaled_expr,"CD86","CD200"))%>%
    dplyr::mutate(CD86_HAVCR2=fn_pairs_score(scaled_expr,"CD86","HAVCR2"))%>%
    dplyr::mutate(TNFRSF14_CD86=fn_pairs_score(scaled_expr,"TNFRSF14","CD86"))%>%
    dplyr::mutate(CTLA4_TNFSF4=fn_pairs_score(scaled_expr,"CTLA4","TNFSF4"))->Project_score
  Project_score
}


all_genes=c("PTEN","CD274","IFNGR1","APLNR","SLC25A13","CPS1","ASL","OTC","SLC25A15","ASS1","IFNG","STAT1","CXCL9","CXCL10","IDO1","HLA-DRA",
            "CD3D","CIITA","CD3E","CCL5","GZMK","CD2","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","GZMB","HLA-A","HLA-F","B2M","TAP1",
            "TAP2","ICAM1","CLECL1","LILRA1","JAK2","CD8A","CD8B","GZMA","PRF1","DKK2","TGFB1","PBRM1","EZH2","KDM1A",
            "VSIR","CD40","PDCD1","CD28","CD80","TNFSF9","CD86","CD276","TNFSF4", "CD27","CD200","HAVCR2","TNFRSF14","CTLA4")
#VSIR == C10orf54
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




#SRP070710-----------------------------------------------------------------------------------------------------------
readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "SRA")%>%
  dplyr::filter(SRA_Study == "SRP070710")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::filter(Response != "NE")%>%
  dplyr::select(Run,Response)->metadata

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_FPKM_expr.txt",header = T)%>%
  dplyr::select(metadata$Run)->SRP070710_expr

rm_ID <- numeric()
for (i in 1:nrow(SRP070710_expr)) {
  if(length(which(is.na(SRP070710_expr[i,]))) >= ncol(SRP070710_expr)*3/4){
    rm_ID <- c(rm_ID,i)
  }
}
SRP070710_expr <- SRP070710_expr[-rm_ID,]

for (i in 1:nrow(SRP070710_expr)) {
  SRP070710_expr[i,which(is.na(SRP070710_expr[i,]))] <- SRP070710_expr[i,]%>%as.numeric()%>%mean(na.rm=TRUE)
}

scale(log2(SRP070710_expr +1),center = TRUE,scale = TRUE)%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% all_genes)%>%
  tibble::column_to_rownames()%>%
  t()%>%
  as.data.frame(stringsAsFactors=FALSE)->scaled_SRP070710_expr

fn_arrange_score_matrix(scaled_SRP070710_expr) ->SRP070710_score


metadata$Response %>%
  gsub("PD","NR",.)%>% gsub("SD","NR",.)%>%
  gsub("PR","R",.)%>%gsub("CR","R",.)->metadata$Response
merge(SRP070710_score,metadata,by.x="rowname",by.y="Run")%>%
  tibble::column_to_rownames()->SRP070710_score


pdf("/data/liull/immune-checkpoint-blockade/pathway_features/SRP070710/SRP070710.pdf")

for (i in 1:(ncol(SRP070710_score)-1)) {
  dplyr::select(SRP070710_score,Response,colnames(SRP070710_score)[i])->single_feature_expr
  
  dplyr::filter(single_feature_expr,Response == "R") %>% dplyr::select(-Response)->R_expr
  dplyr::filter(single_feature_expr,Response == "NR") %>% dplyr::select(-Response)->NR_expr
  wilcox.test(R_expr[,1],NR_expr[,1])$p.value -> anno
  
  if(anno < 0.05){
    print(colnames(SRP070710_score)[i])
    print(i)
    
    #pdf(paste("/data/liull/immune-checkpoint-blockade/pathway_features/","SRP070710/",colnames(SRP070710_score)[i],".pdf",sep = ""))
    ggplot(single_feature_expr)+
      geom_boxplot(aes(x=Response,y=as.numeric(single_feature_expr[,2]),fill=Response),width=0.5)+
      labs(y=colnames(SRP070710_score)[i])+
      annotate(geom="text", x=1.5, y=(max(as.numeric(single_feature_expr[,2]))*3)/4,size=5, label=paste("P.value",signif(anno,2),sep = ":"))+
      scale_fill_manual(values=c("#ef8a62","#67a9cf"))+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border=element_rect(fill = NA),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12))->a
    print(a)
    #grid.newpage()
    
  }
  
}
dev.off()

# [1] "WNT_beta"
# [1] 10
# [1] "TGF_beta"
# [1] 12



#SRP094781-----------------------------------------------------------------------------------------------------------
readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "SRA")%>%
  dplyr::filter(SRA_Study == "SRP094781")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::filter(Response != "NE")%>%
  dplyr::select(Run,Response)->metadata

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_FPKM_expr.txt",header = T)%>%
  dplyr::select(metadata$Run)->SRP094781_expr

rm_ID <- numeric()
for (i in 1:nrow(SRP094781_expr)) {
  if(length(which(is.na(SRP094781_expr[i,]))) >= ncol(SRP094781_expr)*3/4){
    rm_ID <- c(rm_ID,i)
  }
}
SRP094781_expr <- SRP094781_expr[-rm_ID,]

for (i in 1:nrow(SRP094781_expr)) {
  SRP094781_expr[i,which(is.na(SRP094781_expr[i,]))] <- SRP094781_expr[i,]%>%as.numeric()%>%mean(na.rm=TRUE)
}

scale(log2(SRP094781_expr +1),center = TRUE,scale = TRUE)%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% all_genes)%>%
  tibble::column_to_rownames()%>%
  t()%>%
  as.data.frame(stringsAsFactors=FALSE)->scaled_SRP094781_expr

fn_arrange_score_matrix(scaled_SRP094781_expr) ->SRP094781_score


metadata$Response %>%
  gsub("PD","NR",.)%>% gsub("SD","NR",.)%>%
  gsub("PR","R",.)%>%gsub("CR","R",.)->metadata$Response
merge(SRP094781_score,metadata,by.x="rowname",by.y="Run")%>%
  tibble::column_to_rownames()->SRP094781_score

pdf("/data/liull/immune-checkpoint-blockade/pathway_features/SRP094781/SRP094781.pdf")
for (i in 1:(ncol(SRP094781_score)-1)) {
  dplyr::select(SRP094781_score,Response,colnames(SRP094781_score)[i])->single_feature_expr
  
  dplyr::filter(single_feature_expr,Response == "R") %>% dplyr::select(-Response)->R_expr
  dplyr::filter(single_feature_expr,Response == "NR") %>% dplyr::select(-Response)->NR_expr
  wilcox.test(R_expr[,1],NR_expr[,1])$p.value -> anno
  
  if(anno < 0.05){
    print(colnames(SRP094781_score)[i])
    print(i)
    
    #pdf(paste("/data/liull/immune-checkpoint-blockade/pathway_features/","SRP094781/",colnames(SRP094781_score)[i],".pdf",sep = ""))
    ggplot(single_feature_expr)+
      geom_boxplot(aes(x=Response,y=as.numeric(single_feature_expr[,2]),fill=Response),width=0.5)+
      labs(y=colnames(SRP094781_score)[i])+
      annotate(geom="text", x=1.5, y=(max(as.numeric(single_feature_expr[,2]))*3)/4,size=5, label=paste("P.value",signif(anno,2),sep = ":"))+
      scale_fill_manual(values=c("#ef8a62","#67a9cf"))+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border=element_rect(fill = NA),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12))->a
    print(a)
  }
  
}

dev.off()

#/
# [1] "Urea_cycle"
# [1] 4
# [1] "IDO1"
# [1] 11
# [1] "CD40_CD80"
# [1] 19


#SRP150548-----------------------------------------------------------------------------------------------------------
readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "SRA")%>%
  dplyr::filter(SRA_Study == "SRP150548")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::filter(Response != "NE")%>%
  dplyr::filter(Anti_target == "anti-PD1")%>%
  dplyr::select(Run,Response)->metadata

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_FPKM_expr.txt",header = T)%>%
  dplyr::select(metadata$Run)->SRP150548_expr

rm_ID <- numeric()
for (i in 1:nrow(SRP150548_expr)) {
  if(length(which(is.na(SRP150548_expr[i,]))) >= ncol(SRP070710_expr)*3/4){
    rm_ID <- c(rm_ID,i)
  }
}
SRP150548_expr <- SRP150548_expr[-rm_ID,]

for (i in 1:nrow(SRP150548_expr)) {
  SRP150548_expr[i,which(is.na(SRP150548_expr[i,]))] <- SRP150548_expr[i,]%>%as.numeric()%>%mean(na.rm=TRUE)
}

scale(log2(SRP150548_expr +1),center = TRUE,scale = TRUE)%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% all_genes)%>%
  tibble::column_to_rownames()%>%
  t()%>%
  as.data.frame(stringsAsFactors=FALSE)->scaled_SRP150548_expr

fn_arrange_score_matrix(scaled_SRP150548_expr) ->SRP150548_score


metadata$Response %>%
  gsub("PD","NR",.)%>% gsub("SD","NR",.)%>%
  gsub("PR","R",.)%>%gsub("CR","R",.)->metadata$Response
merge(SRP150548_score,metadata,by.x="rowname",by.y="Run")%>%
  tibble::column_to_rownames()->SRP150548_score

for (i in 1:(ncol(SRP150548_score)-1)) {
  
  dplyr::select(SRP150548_score,Response,colnames(SRP150548_score)[i])->single_feature_expr
  if(length(unique(single_feature_expr[,2]))>1){
    dplyr::filter(single_feature_expr,Response == "R") %>% dplyr::select(-Response)->R_expr
    dplyr::filter(single_feature_expr,Response == "NR") %>% dplyr::select(-Response)->NR_expr
    wilcox.test(R_expr[,1],NR_expr[,1])$p.value -> anno
    
    if(anno < 0.05){
      print(colnames(SRP150548_score)[i])
      print(i)
    }else{}
  }
  
}


#gastric cancer ERP107734-----------------------------------------------------------------------------------------------------------
readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "SRA")%>%
  dplyr::filter(SRA_Study == "ERP107734")%>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::filter(Anti_target == "anti-PD1")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::filter(Response != "NE")%>%
  dplyr::select(Run,Response)->metadata

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer_PD1_pretreatment_Symbol_FPKM_expr.txt",header = T)%>%
  dplyr::select(metadata$Run)->ERP107734_expr

rm_ID <- numeric()
for (i in 1:nrow(ERP107734_expr)) {
  if(length(which(is.na(ERP107734_expr[i,]))) >= ncol(ERP107734_expr)*3/4){
    rm_ID <- c(rm_ID,i)
  }
}
ERP107734_expr <- ERP107734_expr[-rm_ID,]

for (i in 1:nrow(ERP107734_expr)) {
  ERP107734_expr[i,which(is.na(ERP107734_expr[i,]))] <- ERP107734_expr[i,]%>%as.numeric()%>%mean(na.rm=TRUE)
}

scale(log2(ERP107734_expr +1),center = TRUE,scale = TRUE)%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% all_genes)%>%
  tibble::column_to_rownames()%>%
  t()%>%
  as.data.frame(stringsAsFactors=FALSE)->scaled_ERP107734_expr

fn_arrange_score_matrix(scaled_ERP107734_expr) ->ERP107734_score


metadata$Response %>%
  gsub("PD","NR",.)%>% gsub("SD","NR",.)%>%
  gsub("PR","R",.)%>%gsub("CR","R",.)->metadata$Response
merge(ERP107734_score,metadata,by.x="rowname",by.y="Run")%>%
  tibble::column_to_rownames()->ERP107734_score

pdf("/data/liull/immune-checkpoint-blockade/pathway_features/ERP107734/ERP107734.pdf")
for (i in 1:(ncol(ERP107734_score)-1)) {
  dplyr::select(ERP107734_score,Response,colnames(ERP107734_score)[i])->single_feature_expr
  
  if(length(unique(single_feature_expr[,2]))>1){
    
    dplyr::filter(single_feature_expr,Response == "R") %>% dplyr::select(-Response)->R_expr
    dplyr::filter(single_feature_expr,Response == "NR") %>% dplyr::select(-Response)->NR_expr
    wilcox.test(R_expr[,1],NR_expr[,1])$p.value -> anno
    
    if(anno < 0.05){
      print(colnames(ERP107734_score)[i])
      print(i)
      
      #pdf(paste("/data/liull/immune-checkpoint-blockade/pathway_features/","ERP107734/",colnames(ERP107734_score)[i],".pdf",sep = ""))
      ggplot(single_feature_expr)+
        geom_boxplot(aes(x=Response,y=as.numeric(single_feature_expr[,2]),fill=Response),width =0.5)+
        labs(y=colnames(ERP107734_score)[i])+
        annotate(geom="text", x=1.5, y=(max(as.numeric(single_feature_expr[,2]))*3)/4,size=5, label=paste("P.value",signif(anno,2),sep = ":"))+
        scale_fill_manual(values=c("#ef8a62","#67a9cf"))+
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border=element_rect(fill = NA),
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12))->a
      print(a)
    }
  }
  
}

dev.off()


#melanoma CTLA4 SRP011540-----------------------------------------------------------------------------------------------------------
readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "dbGAP")%>%
  dplyr::filter(SRA_Study == "SRP011540")%>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::filter(Anti_target == "anti-CTLA4")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::filter(Response != "NE")%>%
  dplyr::filter(Run != "SRR3083584")%>%
  dplyr::select(Run,Second_Response_standard)->metadata
colnames(metadata)[2] <- "Response"

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4_pretreatment_Symbol_FPKM_expr.txt",header = T)%>%
  dplyr::select(metadata$Run)->SRP011540_expr

rm_ID <- numeric()
for (i in 1:nrow(SRP011540_expr)) {
  if(length(which(is.na(SRP011540_expr[i,]))) >= ncol(SRP011540_expr)*3/4){
    rm_ID <- c(rm_ID,i)
  }
}
SRP011540_expr <- SRP011540_expr[-rm_ID,]

for (i in 1:nrow(SRP011540_expr)) {
  SRP011540_expr[i,which(is.na(SRP011540_expr[i,]))] <- SRP011540_expr[i,]%>%as.numeric()%>%mean(na.rm=TRUE)
}

scale(log2(SRP011540_expr +1),center = TRUE,scale = TRUE)%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% all_genes)%>%
  tibble::column_to_rownames()%>%
  t()%>%
  as.data.frame(stringsAsFactors=FALSE)->scaled_SRP011540_expr

fn_arrange_score_matrix(scaled_SRP011540_expr) ->SRP011540_score

metadata$Response %>%
  gsub("long-survival","R",.)->metadata$Response
merge(SRP011540_score,metadata,by.x="rowname",by.y="Run")%>%
  tibble::column_to_rownames()->SRP011540_score

pdf("/data/liull/immune-checkpoint-blockade/pathway_features/SRP011540/SRP011540.pdf")
for (i in 1:(ncol(SRP011540_score)-1)) {
  dplyr::select(SRP011540_score,Response,colnames(SRP011540_score)[i])->single_feature_expr
  
  if(length(unique(single_feature_expr[,2]))>1){
    
    dplyr::filter(single_feature_expr,Response == "R") %>% dplyr::select(-Response)->R_expr
    dplyr::filter(single_feature_expr,Response == "NR") %>% dplyr::select(-Response)->NR_expr
    wilcox.test(R_expr[,1],NR_expr[,1])$p.value -> anno
    
    if(anno < 0.05){
      print(colnames(SRP011540_score)[i])
      print(i)
      
      #pdf(paste("/data/liull/immune-checkpoint-blockade/pathway_features/","SRP011540/",colnames(SRP011540_score)[i],".pdf",sep = ""))
      ggplot(single_feature_expr)+
        geom_boxplot(aes(x=Response,y=as.numeric(single_feature_expr[,2]),fill=Response),width=0.5)+
        labs(y=colnames(SRP011540_score)[i])+
        annotate(geom="text", x=1.5, y=(max(as.numeric(single_feature_expr[,2]))*3)/4,size=5, label=paste("P.value",signif(anno,2),sep = ":"))+
        scale_fill_manual(values=c("#ef8a62","#67a9cf"))+
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border=element_rect(fill = NA),
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12))->a
      print(a)
    }
  }
  
}

dev.off()
