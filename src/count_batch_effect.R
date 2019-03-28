#ComBat--------------------------------------------------------------------------------------
library(magrittr)
library(readr)
library(readxl)
library(dplyr)


#filter melanoma RNA-seq anti-PD1
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(SRA_Study,Run,Response,Biopsy_Time) ->metadata


#expression prepare for batch effect
read.table("/data/liull/immune-checkpoint-blockade/expression/all_count_expression_2.txt",sep="\t",header = T,as.is = TRUE) ->data1
Project=unique(metadata$SRA_Study)
dplyr::filter(metadata,SRA_Study==Project[1]) %>%
  dplyr::select(Run)  %>%
  as.matrix() %>%
  as.character()->Project1_id

dplyr::filter(metadata,SRA_Study==Project[2]) %>%
  dplyr::select(Run)  %>%
  as.matrix() %>%
  as.character()->Project2_id

dplyr::filter(metadata,SRA_Study==Project[3]) %>%
  dplyr::select(Run)  %>%
  as.matrix() %>%
  as.character()->Project3_id

expression1=dplyr::select(data1,gene_id,Project1_id,Project2_id,Project3_id)
#order the expression profile by project to remove batch effect

row.names(expression1)=expression1[,1]
expression1=expression1[,-1]
#make rownames to avoid of sum wrong

Sum_zero=apply(expression1,1,function(x) sum(as.numeric(x),na.rm = TRUE))
IDs_zero=which(Sum_zero==0)
expression2=expression1[-IDs_zero,]#delete the gene has all 0.000 depression

# a=1:nrow(expression3)
# for(j in 1:length(expression3)){
#   which(expression3[,j]<10) ->b
#   intersect(a,b)->a
# }
# expression4=expression3[-a,]#delete the gene has all less than 10 exression


DGEList_expr <- DGEList(counts=expression2)
normalized_expr <- calcNormFactors(DGEList_expr, method="upperquartile")
normalized_loggedCPM_expr = cpm(normalized_expr, log=TRUE, prior.count=2)

#remove batch effect by ComBat
batch1=rep(1,length(Project1_id))
batch2=rep(2,length(Project2_id))
batch3=rep(3,length(Project3_id))
batch=c(batch1,batch2,batch3)
metadata$Response%>%
  gsub("^PD$", "NR",. )%>%
  gsub("^SD$", "NR", .)%>%
  gsub("^PR$", "R", .)%>%
  gsub("^CR$", "R", .)%>%
  gsub("^PRCR$", "R", .)->my_mod
my_mod = model.matrix(~as.factor(my_mod))
combat_edata = ComBat(dat=normalized_loggedCPM_expr, batch=batch, mod=my_mod, par.prior=TRUE, prior.plots=FALSE)

#DEG by limma
dplyr::filter(metadata,Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(Response %in% c("CR","PR","PRCR","R"))-> response
dplyr::filter(metadata,Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(Response %in% c("SD","PD","NR")) -> non_response

dplyr::select(as.data.frame(combat_edata),response$Run,non_response$Run)->ordered_combat_edata
m=1:nrow(ordered_combat_edata)
for(j in 1:ncol(ordered_combat_edata)){
   which(ordered_combat_edata[,j]<0) ->n
   intersect(m,n)->m
 }
expression4=expression3[-a,]#delete the gene has all less than 10 exression

group_list <- factor(c(rep("response",nrow(response)), rep("non-response",nrow(non_response))))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(ordered_combat_edata)


fit <- lmFit(ordered_combat_edata, design)
fit <- eBayes(fit, trend=TRUE)
output <- topTable(fit, coef=2,n=Inf)
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05)->sig_genes

#DEG by 
