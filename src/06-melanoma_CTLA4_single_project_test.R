library(magrittr)

#single one procedure----------------------------------------
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") ->SRA
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="dbGAP") ->dbGAP
rbind(SRA,dbGAP)%>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata
metadata %>% dplyr::filter(Run != "SRR3083584") -> metadata# fastq file 16M
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship

#expression prepare for batch effect
read.table("/data/liull/immune-checkpoint-blockade/expression/all_count_expression_2.txt",sep="\t",header = T,as.is = TRUE) ->all_count_data
Project=unique(metadata$SRA_Study)

i=2
dplyr::filter(metadata,SRA_Study==Project[i])%>%
  dplyr::select(Run,Response)->traits
dplyr::filter(traits,Response %in% c("CR","PR","R","X")) ->response
dplyr::filter(traits,Response %in% c("SD","PD","NR")) ->non_response
  
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
  
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC>1)->up
tibble::rownames_to_column(output) %>% dplyr::filter(P.Value<0.05) %>% dplyr::filter(logFC< -1)->down
  
merge(relationship,up,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% up$rowname) ->up2
up2[order(up2$logFC,decreasing = TRUE),]->up2
merge(relationship,down,by.x="Ensembl_ID",by.y="rowname",all=TRUE)%>%
  dplyr::filter(Ensembl_ID %in% down$rowname) ->down2
down2[order(down2$logFC),]->down2
dim(up2)
dim(down2)
write.table(up2,paste("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/single_test/",Project[i],"_all_up.txt",sep="" ),quote = FALSE,row.names = FALSE,col.names = TRUE)
write.table(down2,paste("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/single_test/",Project[i],"_all_down.txt",sep="" ),quote = FALSE,row.names = FALSE,col.names = TRUE)

#single result compare with all_together-------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/single_test/SRP011540_all_up.txt",header = T,as.is = TRUE) ->SRP011540_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/single_test/SRP011540_all_down.txt",header = T,as.is = TRUE) ->SRP011540_down

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_all_up.txt",header = T,as.is = TRUE) ->CTLA4_all_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_all_down.txt",header = T,as.is = TRUE) ->CTLA4_all_down  

venn.diagram(list(SRP011540=SRP011540_up$Ensembl_ID,CTLA4_all=CTLA4_all_up$rowname),filename=NULL,fill=c("cornflowerblue", "yellow"))->venn_up
grid.draw(venn_up)
ggsave(
  filename = 'venn_up_interaction.pdf',
  plot = venn_up,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/single_test/',
  width = 12,
  height = 8
)

venn.diagram(list(SRP011540=SRP011540_down$Ensembl_ID,CTLA4_all=CTLA4_all_down$rowname),filename=NULL,fill=c("cornflowerblue", "yellow"))->venn_down
grid.draw(venn_down)
ggsave(
  filename = 'venn_down_interaction.pdf',
  plot = venn_down,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/single_test/',
  width = 12,
  height = 8
)







