#use the all enriched functions' gene as features(GO,KEGG)
library()


#GO interaction-------------------------------------------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/up_enrichGO.txt",
           header = T,as.is = TRUE,sep = "\t",quote = "")%>%
  dplyr::filter(p.adjust <= 0.01)->GO_melanoma_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/up_enrichGO.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->GO_melanoma_CTLA4_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_enrichGO.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->GO_gastric_PD1_up
intersect(GO_melanoma_PD1_up$Description,GO_melanoma_CTLA4_up$Description) %>% intersect(GO_gastric_PD1_up$Description)->up_Intersection

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/down_enrichGO.txt",
           header = T,as.is = TRUE,sep = "\t",quote = "")%>%
  dplyr::filter(p.adjust <= 0.01)->GO_melanoma_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/down_enrichGO.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->GO_melanoma_CTLA4_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_enrichGO.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->GO_gastric_PD1_down
intersect(GO_melanoma_PD1_down$Description,GO_melanoma_CTLA4_down$Description) %>% intersect(GO_gastric_PD1_down$Description)->down_Intersection

#KEGG interaction-----------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/up_enrichKEGG.txt",
           header = T,as.is = TRUE,sep = "\t",quote = "")%>%
  dplyr::filter(p.adjust <= 0.01)->KEGG_melanoma_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/up_enrichKEGG.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->KEGG_melanoma_CTLA4_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_enrichKEGG.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->KEGG_gastric_PD1_up
intersect(KEGG_melanoma_PD1_up$Description,KEGG_melanoma_CTLA4_up$Description) %>% intersect(KEGG_gastric_PD1_up$Description)->KEGG_up_Intersection



read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/down_enrichKEGG.txt",
           header = T,as.is = TRUE,sep = "\t",quote = "")%>%
  dplyr::filter(p.adjust <= 0.01)->KEGG_melanoma_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/down_enrichKEGG.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->KEGG_melanoma_CTLA4_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_enrichKEGG.txt",
           header = T,as.is = TRUE,sep = "\t")%>%
  dplyr::filter(p.adjust <= 0.01)->KEGG_gastric_PD1_down
intersect(KEGG_melanoma_PD1_down$Description,KEGG_melanoma_CTLA4_down$Description) %>% intersect(KEGG_gastric_PD1_down$Description)->KEGG_down_Intersection
#0

#get all genes
fn_all_pathway_genes_GO <- function(enrich_result,GO_Description){
  
  dplyr::filter(enrich_result,Description %in% GO_Description)%>%
    dplyr::select(geneID)-> Genes
  
  sapply(Genes$geneID,function(x) strsplit(x,"/"))%>%
    unlist()%>%
    as.character()%>%
    unique()
  
}
readxl::read_excel("/data/liull/reference/All_EntrezID_Symbl_NCBI.xlsx",col_names = TRUE)%>%as.data.frame(stringsAsFactors)->EntrezID_Symbl
fn_all_pathway_genes_KEGG <- function(enrich_result,KEGG_Description){
  
  dplyr::filter(enrich_result,Description %in% KEGG_Description)%>%
    dplyr::select(geneID)-> Genes
  
  sapply(Genes$geneID,function(x) strsplit(x,"/"))%>%
    unlist()%>%
    as.character()%>%
    unique()->Genes
  data.frame(EntrezID=Genes)%>%
    merge(EntrezID_Symbl,by.x="EntrezID",by.y="GeneID")->Genes
  Genes$Symbol
  
}

gene_set=list()
gene_set$T_cell_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"T cell activation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"T cell activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"T cell activation"))
gene_set$regulation_of_T_cell_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of T cell activation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of T cell activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of T cell activation"))
gene_set$positive_regulation_of_T_cell_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of T cell activation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of T cell activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of T cell activation"))
#
#
gene_set$T_cell_costimulation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"T cell costimulation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"T cell costimulation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"T cell costimulation"))
#
#
gene_set$lymphocyte_costimulation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"lymphocyte costimulation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"lymphocyte costimulation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"lymphocyte costimulation"))
gene_set$lymphocyte_differentiation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"lymphocyte differentiation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"lymphocyte differentiation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"lymphocyte differentiation"))
gene_set$regulation_of_lymphocyte_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of lymphocyte activation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of lymphocyte activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of lymphocyte activation"))
gene_set$positive_regulation_of_lymphocyte_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of lymphocyte activation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of lymphocyte activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of lymphocyte activation"))
gene_set$lymphocyte_proliferation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"lymphocyte proliferation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"lymphocyte proliferation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"lymphocyte proliferation"))
gene_set$regulation_of_lymphocyte_proliferation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of lymphocyte proliferation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of lymphocyte proliferation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of lymphocyte proliferation"))
#
#
gene_set$leukocyte_cell_cell_adhesion <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"leukocyte cell-cell adhesion"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"leukocyte cell-cell adhesion"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"leukocyte cell-cell adhesion"))
gene_set$leukocyte_differentiation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"leukocyte differentiation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"leukocyte differentiation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"leukocyte differentiation"))
gene_set$leukocyte_proliferation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"leukocyte proliferation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"leukocyte proliferation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"leukocyte proliferation"))
gene_set$positive_regulation_of_leukocyte_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of leukocyte activation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of leukocyte activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of leukocyte activation"))
gene_set$positive_regulation_of_leukocyte_cell_cell_adhesion <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of leukocyte cell-cell adhesion"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of leukocyte cell-cell adhesion"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of leukocyte cell-cell adhesion"))
gene_set$regulation_of_leukocyte_proliferation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of leukocyte proliferation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of leukocyte proliferation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of leukocyte proliferation"))
gene_set$regulation_of_leukocyte_cell_cell_adhesion <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of leukocyte cell-cell adhesion"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of leukocyte cell-cell adhesion"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of leukocyte cell-cell adhesion"))
#
#
gene_set$interferon_gamma_production <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"interferon-gamma production"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"interferon-gamma production"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"interferon-gamma production"))
gene_set$regulation_of_interferon_gamma_production <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of interferon-gamma production"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of interferon-gamma production"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of interferon-gamma production"))
#
#
gene_set$mononuclear_cell_proliferation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"mononuclear cell proliferation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"mononuclear cell proliferation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"mononuclear cell proliferation"))
gene_set$regulation_of_mononuclear_cell_proliferation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of mononuclear cell proliferation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of mononuclear cell proliferation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of mononuclear cell proliferation"))
#
#
gene_set$MHC_class_II_protein_complex <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"MHC class II protein complex"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"MHC class II protein complex"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"MHC class II protein complex"))
#
#
gene_set$regulation_of_cell_cell_adhesion <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of cell-cell adhesion"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of cell-cell adhesion"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of cell-cell adhesion"))
gene_set$positive_regulation_of_cell_adhesion <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of cell adhesion"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of cell adhesion"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of cell adhesion"))
gene_set$positive_regulation_of_cell_cell_adhesion <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of cell-cell adhesion"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of cell-cell adhesion"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of cell-cell adhesion"))
#
gene_set$antigen_receptor_mediated_signaling_pathway <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"antigen receptor-mediated signaling pathway"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"antigen receptor-mediated signaling pathway"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"antigen receptor-mediated signaling pathway"))
gene_set$cellular_defense_response <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"cellular defense response"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"cellular defense response"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"cellular defense response"))
gene_set$immune_response_activating_cell_surface_receptor_signaling_pathway <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"immune response-activating cell surface receptor signaling pathway"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"immune response-activating cell surface receptor signaling pathway"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"immune response-activating cell surface receptor signaling pathway"))
gene_set$immune_response_regulating_cell_surface_receptor_signaling_pathway <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"immune response-regulating cell surface receptor signaling pathway"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"immune response-regulating cell surface receptor signaling pathway"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"immune response-regulating cell surface receptor signaling pathway"))
gene_set$immunological_synapse <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"immunological synapse"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"immunological synapse"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"immunological synapse"))
gene_set$positive_regulation_of_cell_activation <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"positive regulation of cell activation"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"positive regulation of cell activation"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"positive regulation of cell activation"))
gene_set$regulation_of_immune_effector_process <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"regulation of immune effector process"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"regulation of immune effector process"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"regulation of immune effector process"))
gene_set$response_to_lipopolysaccharide <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"response to lipopolysaccharide"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"response to lipopolysaccharide"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"response to lipopolysaccharide"))
gene_set$response_to_molecule_of_bacterial_origin <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"response to molecule of bacterial origin"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"response to molecule of bacterial origin"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"response to molecule of bacterial origin"))
gene_set$side_of_membrane <- union(fn_all_pathway_genes_GO(GO_melanoma_PD1_up,"side of membrane"),fn_all_pathway_genes_GO(GO_melanoma_CTLA4_up,"side of membrane"))%>%
  union(fn_all_pathway_genes_GO(GO_gastric_PD1_up,"side of membrane"))
#
#
gene_set$Hematopoietic_cell_lineage <- union(fn_all_pathway_genes_KEGG(KEGG_melanoma_PD1_up,"Hematopoietic cell lineage"),fn_all_pathway_genes_KEGG(KEGG_melanoma_CTLA4_up,"Hematopoietic cell lineage"))%>%
  union(fn_all_pathway_genes_KEGG(KEGG_gastric_PD1_up,"Hematopoietic cell lineage"))
gene_set$Staphylococcus_aureus_infection <- union(fn_all_pathway_genes_KEGG(KEGG_melanoma_PD1_up,"Staphylococcus aureus infection"),fn_all_pathway_genes_KEGG(KEGG_melanoma_CTLA4_up,"Staphylococcus aureus infection"))%>%
  union(fn_all_pathway_genes_KEGG(KEGG_gastric_PD1_up,"Staphylococcus aureus infection"))
gene_set$Intestinal_immune_network_for_IgA_production <- union(fn_all_pathway_genes_KEGG(KEGG_melanoma_PD1_up,"Intestinal immune network for IgA production"),fn_all_pathway_genes_KEGG(KEGG_melanoma_CTLA4_up,"Intestinal immune network for IgA production"))%>%
  union(fn_all_pathway_genes_KEGG(KEGG_gastric_PD1_up,"Intestinal immune network for IgA production"))
gene_set$CAMs <- union(fn_all_pathway_genes_KEGG(KEGG_melanoma_PD1_up,"Cell adhesion molecules (CAMs)"),fn_all_pathway_genes_KEGG(KEGG_melanoma_CTLA4_up,"Cell adhesion molecules (CAMs)"))%>%
  union(fn_all_pathway_genes_KEGG(KEGG_gastric_PD1_up,"Cell adhesion molecules (CAMs)"))
gene_set$Th17_cell_differentiation <- union(fn_all_pathway_genes_KEGG(KEGG_melanoma_PD1_up,"Th17 cell differentiation"),fn_all_pathway_genes_KEGG(KEGG_melanoma_CTLA4_up,"Th17 cell differentiation"))%>%
  union(fn_all_pathway_genes_KEGG(KEGG_gastric_PD1_up,"Th17 cell differentiation"))


fn_get_feature_score <- function(expr,gene_list){
  
  score_frame <- matrix(nrow = ncol(expr),ncol = length(gene_list))
  colnames(score_frame) = names(gene_list)
  rownames(score_frame) = colnames(expr)
  
  t(expr)%>%as.data.frame(stringsAsFactors=FALSE)->t_expr
  
  for (i in 1:ncol(score_frame)) {
    
    which(names(gene_list)[]==colnames(score_frame)[i]) -> feature_id
    gene_list[feature_id] %>% unlist()%>% as.character()->feature_gene
    
    t_expr %>% dplyr::select(feature_gene) -> feature_gene_expr
    apply(feature_gene_expr, 1, mean)->feature_score
    names(feature_score) <- NULL
    
    score_frame[,i] <- feature_score
    
  }
  score_frame
}

#train in melanoma PD1------------------------------

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
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_log2CPM_expr.txt",
           header = T)%>%
  dplyr::select(metadata$Run)%>%
  tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% unique(as.character(unlist(gene_set))))%>%
  tibble::column_to_rownames()->melanoma_PD1_log2CPM

# scale(melanoma_PD1_log2CPM,center = TRUE,scale = TRUE)%>%
#   as.data.frame(stringsAsFactors=FALSE)->scaled_melanoma_PD1_log2CPM
fn_get_feature_score(melanoma_PD1_log2CPM,gene_set)%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  tibble::rownames_to_column()%>%
  merge(metadata[,c(1,3)],by.x="rowname",by.y="Run")%>%
  tibble::column_to_rownames()-> train_features_score
train_features_score$Response=as.factor(train_features_score$Response)

#test in melanoma CTLA4------------------------------

readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "dbGAP")%>%
  dplyr::filter(Cancer == "melanoma")%>%
  dplyr::filter(Anti_target =="anti-CTLA4")%>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(Run,SRA_Study,Second_Response_standard)%>%
  dplyr::filter(Run != "SRR3083584")->metadata
metadata$Second_Response_standard %>%
  gsub("long-survival","R",.)->metadata$Second_Response_standard
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4_pretreatment_Symbol_log2CPM_expr.txt",
           header = T)%>%
  dplyr::select(metadata$Run)%>%
  tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% unique(as.character(unlist(gene_set))))%>%
  tibble::column_to_rownames()->melanoma_CTLA4_log2CPM

# scale(melanoma_CTLA4_log2CPM,center = TRUE,scale = TRUE)%>%
#   as.data.frame(stringsAsFactors=FALSE)->scaled_melanoma_CTLA4_log2CPM
fn_get_feature_score(melanoma_CTLA4_log2CPM,gene_set)%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  tibble::rownames_to_column()%>%
  merge(metadata[,c(1,3)],by.x="rowname",by.y="Run")%>%
  tibble::column_to_rownames()-> test_features_score
test_features_score$Second_Response_standard=as.factor(test_features_score$Second_Response_standard)

#randomforest----------
ROC=c()
NR_1000_id=list()
for (i in 1:1000) {
  
  which(train_features_score$Response[]=="NR") %>% sample(26)-> NR_id
  which(train_features_score$Response[]=="R") -> R_id
  
  
  randomForest(x=train_features_score[c(NR_id,R_id),-41], y=train_features_score[c(NR_id,R_id),]$Response, ntree = 500,type="classification",
               importance=T )->features_rf
  
  pred_test=predict(features_rf,test_features_score[,-41])
  test_roc <- roc(test_features_score$Second_Response_standard,as.numeric(pred_test))
  ROC[i] <- test_roc$auc[1]
  NR_1000_id[[i]] <- NR_id
  
}
NR_1000_id[[which(ROC[]==max(ROC))]]->NR_id
#gastric cancer
readxl::read_xlsx("/data/liucj/data/immune-checkpoint-blockade/all_metadata_available.xlsx",sheet = "SRA")%>%
  dplyr::filter(Cancer == "gastric cancer")%>%
  dplyr::filter(Anti_target =="anti-PD1")%>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::filter(Biopsy_Time == "pre-treatment")%>%
  dplyr::select(Run,SRA_Study,Response)->metadata
metadata$Response %>%
  gsub("PD","NR",.)%>% gsub("SD","NR",.)%>%
  gsub("PR","R",.)%>%gsub("CR","R",.)->metadata$Response
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer_PD1_pretreatment_Symbol_log2CPM_expr.txt",
           header = T)%>%
  dplyr::select(metadata$Run)%>%
  tibble::rownames_to_column()%>%
  dplyr::filter(rowname %in% unique(as.character(unlist(gene_set))))%>%
  tibble::column_to_rownames()->gastric_cancer_PD1_log2CPM

# scale(gastric_cancer_PD1_log2CPM,center = TRUE,scale = TRUE)%>%
#   as.data.frame(stringsAsFactors=FALSE)->scaled_gastric_cancer_PD1_log2CPM
fn_get_feature_score(gastric_cancer_PD1_log2CPM,gene_set)%>%
  as.data.frame(stringsAsFactors=FALSE)%>%
  tibble::rownames_to_column()%>%
  merge(metadata[,c(1,3)],by.x="rowname",by.y="Run")%>%
  tibble::column_to_rownames()-> gastric_cancer_features_score
gastric_cancer_features_score$Response=as.factor(gastric_cancer_features_score$Response)

pred_gastric=predict(features_rf,gastric_cancer_features_score[,-41])
test_roc <- roc(gastric_cancer_features_score$Response,as.numeric(pred_gastric))
