library(magrittr)
library(edgeR)
library(VennDiagram)

#single one procedure----------------------------------------
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship

#expression prepare for batch effect
read.table("/data/liull/immune-checkpoint-blockade/expression/all_count_expression_2.txt",sep="\t",header = T,as.is = TRUE) ->all_count_data
Project=unique(metadata$SRA_Study)

#SRP070710_standard vs SRP070710_me--------------------------------------------------------------------------
i=1
dplyr::filter(metadata,SRA_Study==Project[i])%>%
  dplyr::select(Run,Response)->traits
dplyr::filter(traits,Response %in% c("CR","PR","R")) ->response
dplyr::filter(traits,Response %in% c("SD","PD","NR")) ->non_response

all_count_data%>%
  dplyr::select(gene_id,response$Run,non_response$Run)%>%
  dplyr::filter(gene_id %in% relationship$Ensembl_ID) %>%
  merge(relationship,.,by.x="Ensembl_ID",by.y="gene_id")%>%
  dplyr::select(-Ensembl_ID,-GeneID)->single_expr  #symbol-PD1-pretreatment-log2CPM-expr
dim(single_expr)#[1] 22114    28

factors=factor(single_expr$Symbol)
merged_expression=tapply(single_expr[,2],factors,median)
for (i in 3:ncol(single_expr)) {
  temp=tapply(single_expr[,i],factors,median)
  merged_expression=cbind(merged_expression,temp)
}
colnames(merged_expression)=colnames(single_expr)[2:ncol(single_expr)]  #trans ensembl id to symbol and merged
dim(merged_expression)#[1] 22094    27


  
DGEList_expr <- DGEList(counts=merged_expression)
normalized_expr <- calcNormFactors(DGEList_expr, method="upperquartile")
normalized_loggedCPM_expr = cpm(normalized_expr, log=TRUE, prior.count=2)
  
keep <- rowSums(normalized_loggedCPM_expr>0) >= 2
normalized_loggedCPM_expr <- normalized_loggedCPM_expr[keep,]
dim(normalized_loggedCPM_expr)#[1] 15405    27
#delete the gene has less than 2 sample exression CPM<1(log2CPM<0)
  
group_list <- factor(c(rep("response",nrow(response)), rep("non_response",nrow(non_response))))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(normalized_loggedCPM_expr)
  
fit <- lmFit(normalized_loggedCPM_expr, design)
fit2 <- eBayes(fit)
output <- topTable(fit2, coef=2, n=Inf)
  
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.1) %>% dplyr::filter(logFC>log2(1.5))->up#372
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.1) %>% dplyr::filter(logFC<log2(1/1.5))->down#1072

write.table(up,paste("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP070710_all_up.txt",sep="" ),quote = FALSE,row.names = FALSE,col.names = TRUE)
write.table(down,paste("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP070710_all_down.txt",sep="" ),quote = FALSE,row.names = FALSE,col.names = TRUE)


#up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP070710_standard_up.txt",sep="\t",header = T,as.is = TRUE)%>%
  as.matrix()%>%
  as.character()->SRP070710_standard_up

venn.diagram(list(SRP070710=up$rowname, SRP070710_standard=SRP070710_standard_up),main="SRP070710 VS SRP070710_standard: up",
             cex=2,margin = 0.1,imagetype = "png",
             filename = "/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP070710_venn_up.png",
             fill=c("red","yellow"),na="remove")

#down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP070710_standard_down.txt",sep="\t",header = T,as.is = TRUE)%>%
  as.matrix()%>%
  as.character()->SRP070710_standard_down

venn.diagram(list(SRP070710=down$rowname, SRP070710_standard=SRP070710_standard_down),main="SRP070710 VS SRP070710_standard: down",
             cex=2,margin = 0.1,imagetype = "png",
             filename = "/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP070710_venn_down.png",
             fill=c("red","yellow"),na="remove")


#SRP094781_standard vs SRP094781_me -------------------------------------------
i=3
dplyr::filter(metadata,SRA_Study==Project[i])%>%
  dplyr::select(Run,Response)->traits
dplyr::filter(traits,Response %in% c("CR","PR")) ->response
dplyr::filter(traits,Response %in% c("PD","SD")) ->non_response

dplyr::select(all_count_data,gene_id,response$Run,non_response$Run)->single_expr
row.names(single_expr)=single_expr[,1]
single_expr=single_expr[,-1]

DGEList_expr <- DGEList(counts=single_expr)
normalized_expr <- calcNormFactors(DGEList_expr, method="upperquartile")
normalized_loggedCPM_expr = cpm(normalized_expr, log=TRUE, prior.count=2)

keep <- rowSums(normalized_loggedCPM_expr>0) >= 2
normalized_loggedCPM_expr <- normalized_loggedCPM_expr[keep,]
dim(normalized_loggedCPM_expr)
#delete the gene has less than 2 sample exression CPM<1(log2CPM<0)

group_list <- factor(c(rep("response",nrow(response)), rep("non_response",nrow(non_response))))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(normalized_loggedCPM_expr)

fit <- lmFit(normalized_loggedCPM_expr, design)
fit2 <- eBayes(fit)
output <- topTable(fit2, coef=2, n=Inf)

tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.1) %>% dplyr::filter(logFC>log2(1.2))->up
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.1) %>% dplyr::filter(logFC<log2(1/1.2))->down

merge(relationship,up,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% up$rowname) ->up2
up2[order(up2$logFC,decreasing = TRUE),]->up2
merge(relationship,down,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% down$rowname) ->down2
down2[order(down2$logFC),]->down2
dim(up2)
dim(down2)
write.table(up2,paste("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/",Project[i],"_all_up_1.2.txt",sep="" ),quote = FALSE,row.names = FALSE,col.names = TRUE)
write.table(down2,paste("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/",Project[i],"_all_down_1.2.txt",sep="" ),quote = FALSE,row.names = FALSE,col.names = TRUE)

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP094781_all_up_1.2.txt",header = T,as.is = TRUE) ->SRP094781_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP094781_all_down_1.2.txt",header = T,as.is = TRUE) ->SRP094781_down
#up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP094781_standard_up.txt",sep="\t",header = T,as.is = TRUE)%>%
  as.matrix()%>%
  as.character()->SRP094781_standard_up

dplyr::filter(SRP094781_up,Ensembl_ID %in% grep("ENSG",SRP094781_up$Ensembl_ID,value=T))%>%
  dplyr::select(Symbol)%>%
  as.matrix()%>%
  as.character()->SRP094781_me_up

setwd("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/")
venn.diagram(list(SRP094781=SRP094781_me_up, SRP094781_standard=SRP094781_standard_up),main="SRP094781 VS SRP094781_standard: up",cex=2,margin = 0.1,imagetype = "png",filename="SRP094781_venn_up.png",fill=c("green","yellow"),na="remove")


#down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP094781_standard_down.txt",sep="\t",header = T,as.is = TRUE)%>%
  as.matrix()%>%
  as.character()->SRP094781_standard_down

dplyr::filter(SRP094781_down,Ensembl_ID %in% grep("ENSG",SRP094781_down$Ensembl_ID,value=T))%>%
  dplyr::select(Symbol)%>%
  as.matrix()%>%
  as.character()->SRP094781_me_down

venn.diagram(list(SRP094781=SRP094781_me_down, SRP094781_standard=SRP094781_standard_down),main="SRP094781 VS SRP094781_standard: down",cex=2,margin = 0.1,imagetype = "png",filename="SRP094781_venn_down.png",fill=c("green","yellow"),na="remove")


#SRP150548's paper no DEG--------------




#single result compare with all_together-------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP070710_all_up_1.5.txt",header = T,as.is = TRUE) ->SRP070710_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP070710_all_down_1.5.txt",header = T,as.is = TRUE) ->SRP070710_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP094781_all_up_1.5.txt",header = T,as.is = TRUE) ->SRP094781_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP094781_all_down_1.5.txt",header = T,as.is = TRUE) ->SRP094781_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP150548_all_up_1.5.txt",header = T,as.is = TRUE) ->SRP150548_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/SRP150548_all_down_1.5.txt",header = T,as.is = TRUE) ->SRP150548_down

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_all_up.txt",header = T,as.is = TRUE) ->PD1_all_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_all_down.txt",header = T,as.is = TRUE) ->PD1_all_down

venn.diagram(list(SRP070710=SRP070710_up$Ensembl_ID, SRP094781=SRP094781_up$Ensembl_ID,SRP150548=SRP150548_up$Ensembl_ID,PD1_all=PD1_all_up$rowname),filename=NULL,fill=c("cornflowerblue", "green", "yellow","red"))->venn_up
grid.draw(venn_up)
ggsave(
  filename = 'venn_up_interaction.pdf',
  plot = venn_up,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/',
  width = 12,
  height = 8
)

# > length(intersect(SRP150548_up$Ensembl_ID,PD1_all_up$rowname))/length(SRP150548_up$Ensembl_ID)
# [1] 0.03532609
# > length(intersect(SRP070710_up$Ensembl_ID,PD1_all_up$rowname))/length(SRP070710_up$Ensembl_ID)
# [1] 0.1541667
# > length(intersect(SRP094781_up$Ensembl_ID,PD1_all_up$rowname))/length(SRP094781_up$Ensembl_ID)
# [1] 0.2745098




venn.diagram(list(SRP070710=SRP070710_down$Ensembl_ID, SRP094781=SRP094781_down$Ensembl_ID,SRP150548=SRP150548_down$Ensembl_ID,PD1_all=PD1_all_down$rowname),filename=NULL,fill=c("cornflowerblue", "green", "yellow","red"))->venn_down
grid.draw(venn_down)
ggsave(
  filename = 'venn_down_interaction.pdf',
  plot = venn_down,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/single_test/',
  width = 12,
  height = 8
)

# > length(intersect(SRP150548_down$Ensembl_ID,PD1_all_down$rowname))/length(SRP150548_down$Ensembl_ID)
# [1] 0.0195713
# > length(intersect(SRP070710_down$Ensembl_ID,PD1_all_down$rowname))/length(SRP070710_down$Ensembl_ID)
# [1] 0.1036364
# > length(intersect(SRP094781_down$Ensembl_ID,PD1_all_down$rowname))/length(SRP094781_down$Ensembl_ID)
# [1] 0.3048128

