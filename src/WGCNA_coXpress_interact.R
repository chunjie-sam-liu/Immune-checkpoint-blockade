library(magrittr)
#WGCNA-----------------------------------------------------------------------------------------------------------------------------------
read.table("/data/liull/test_WGCNA_R/module.color.txt",header = F,as.is = TRUE,skip = 1) ->R_color_module #response modules
cbind(R_color_module[1:2,],R_color_module[3:4,])%>%
  cbind(R_color_module[5:6,])%>%
  cbind(R_color_module[7:8,])%>%
  cbind(R_color_module[9:10,])%>%
  cbind(R_color_module[11:12,])%>%
  cbind(R_color_module[13:14,])%>%
  cbind(R_color_module[15:16,])%>%
  t()%>%
  as.data.frame(stringsAsFactors=FALSE)->R_color_module
colnames(R_color_module)=c("module","Num")
rownames(R_color_module)=NULL
dplyr::filter(R_color_module,module != "grey")->R_color_module

read.table("/data/liull/test_WGCNA_R/raw_module.assign.txt",header = T,as.is = TRUE) %>%
  dplyr::filter(module != "grey")->R_color_gene

R_list_sets=list()
R_module_names=character()
for (i in 1:nrow(R_color_module)) {
  
  R_color_module$module[i] -> R_module_names[i]
  dplyr::filter(R_color_gene,module == R_module_names[i])%>%
    dplyr::select(gene)%>%
    as.matrix()%>%
    as.character()%>%
    list()->R_list_sets[i]
  
}
names(R_list_sets)=paste("R_",R_module_names,sep="")


read.table("/data/liull/test_WGCNA_NR/module.color.txt",header = F,as.is = TRUE,skip = 1,fill = T) ->NR_color_module #non_response modules
cbind(NR_color_module[1:2,],NR_color_module[3:4,])%>%
  cbind(NR_color_module[5:6,])%>%
  cbind(NR_color_module[7:8,])%>%
  cbind(NR_color_module[9:10,])%>%
  cbind(NR_color_module[11:12,])%>%
  cbind(NR_color_module[13:14,])%>%
  cbind(NR_color_module[15:16,])%>%
  cbind(NR_color_module[17:18,])%>%
  cbind(NR_color_module[19:20,])%>%
  t()%>%
  as.data.frame(stringsAsFactors=FALSE)->NR_color_module
colnames(NR_color_module)=c("module","Num")
rownames(NR_color_module)=NULL
dplyr::filter(NR_color_module,module != "grey") %>%
  dplyr::filter(module != "")->NR_color_module

read.table("/data/liull/test_WGCNA_NR/raw_module.assign.txt",header = T,as.is = TRUE) %>%
  dplyr::filter(module != "grey") ->NR_color_gene

NR_list_sets=list()
NR_module_names=character()
for (i in 1:nrow(NR_color_module)) {
  
  NR_color_module$module[i] -> NR_module_names[i]
  dplyr::filter(NR_color_gene,module == NR_module_names[i])%>%
    dplyr::select(gene)%>%
    as.matrix()%>%
    as.character()%>%
    list()->NR_list_sets[i]
  
}
names(NR_list_sets)=paste("NR_",NR_module_names,sep="")


list_sets_WGCNA=c(R_list_sets,NR_list_sets)

#coxpress------------------------------------------------------------------------------------------------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/pre_PD1_filtered_symbol_expr.txt",header = T,as.is = TRUE) ->pre_PD1_expr
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response,Biopsy_Time,Survival_time,Survival_status,Age,Gender) ->metadata
dplyr::filter(metadata,Response %in% c("CR","PR","R"))-> response#26
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response#59
dplyr::select(pre_PD1_expr,response$Run,non_response$Run)->ordered_PD1_expr
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/selected_groups.txt",header = T,as.is = TRUE)->selected_groups

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/g.txt",header = F,as.is = TRUE) ->g_dataframe
g=g_dataframe$V2
names(g)=g_dataframe$V1

list_sets_coXpress=list()
module_names=character()
for (i in 1:nrow(selected_groups)) {
  
  module=selected_groups$group[i]
  inspect.group(ordered_PD1_expr,g,selected_groups$group[i],1:26, 27:85)->sets
  list_sets_coXpress[i]=list(union(sets$GeneA,sets$GeneB))
  module_names[i]=paste("R_module",module,sep="_")
  
}
names(list_sets_coXpress)=module_names


for (i in 1:length(list_sets_WGCNA)) {
  for (j in 1:length(list_sets_coXpress)) {
    Interaction=intersect(list_sets_WGCNA[[i]],list_sets_coXpress[[j]])
    # length(intersect(list_sets_WGCNA[[i]],list_sets_coXpress[[j]]))/length(list_sets_WGCNA[[i]])->ratio1
    # length(intersect(list_sets_WGCNA[[i]],list_sets_coXpress[[j]]))/length(list_sets_coXpress[[j]])->ratio2
    if(length(Interaction)!=0){
      print(Interaction)
    
    }
     
  }
}


# [1] 35
# [1] 4
# 1 interaction:"ABCG5"

