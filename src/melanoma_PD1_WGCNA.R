#WGCNA for melanoma PD1
library(WGCNA)
library(magrittr)
enableWGCNAThreads(nThreads =20 )
#prepared data traits----------------------------------------------------------------------------
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::filter(Response != "NE")%>%
  dplyr::select(Run,Response)%>%
  as.data.frame()->melanoma_PD1
melanoma_PD1$Response%>%
  gsub("^PD$", "NR",. )%>%
  gsub("^SD$", "NR", .)%>%
  gsub("^PR$", "R", .)%>%
  gsub("^CR$", "R", .)->melanoma_PD1$Response#85 samples
datTraits=as.data.frame(melanoma_PD1$Response)
rownames(datTraits)=melanoma_PD1$Run
colnames(datTraits)="Response"   


#prepared data expr------------------------------------------------------------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/PD1_removed_batch_expression.txt",header = T,as.is = TRUE) ->all_expression
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",header = T,as.is = TRUE,sep="\t") -> relationship

# tibble::rownames_to_column(all_expression) %>%
#   dplyr::select(rowname,melanoma_PD1$Run)->m_PD1_pre
# keep <- rowSums(m_PD1_pre>0) >= 2
# m_PD1_pre <- m_PD1_pre[keep,]
# write.table(m_PD1_pre,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/melanoma_PD1_pretreatment_ENS_NONH_log2CPM_largerthan_2_expr.txt",row.names = TRUE,col.names = TRUE,quote=FALSE,sep="\t")

tibble::rownames_to_column(all_expression) %>%
  dplyr::select(rowname,melanoma_PD1$Run)%>%
  dplyr::filter(rowname %in% relationship$Ensembl_ID) %>%
  merge(relationship,.,by.x="Ensembl_ID",by.y="rowname")%>%
  dplyr::select(-Ensembl_ID,-GeneID)->melanoma_PD1_expr0  #symbol-PD1-pretreatment-log2CPM-expr

factors=factor(melanoma_PD1_expr0$Symbol)
merged_expression=tapply(melanoma_PD1_expr0[,2],factors,median)
for (i in 3:ncol(melanoma_PD1_expr0)) {
  temp=tapply(melanoma_PD1_expr0[,i],factors,median)
  merged_expression=cbind(merged_expression,temp)
}
colnames(merged_expression)=colnames(melanoma_PD1_expr0)[2:ncol(melanoma_PD1_expr0)]  #trans ensembl id to symbol and merged

keep <- rowSums(merged_expression>0) >= 2
melanoma_PD1_expr1 <- merged_expression[keep,]  #keep the gene has more than 2 CPM>2's sample   [1] 16378    85
write.table(melanoma_PD1_expr1,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/melanoma_PD1_pretreatment_symbol_log2CPM_largerthan_2_expr.txt",row.names = TRUE,col.names = TRUE,quote=FALSE,sep="\t")

expr = as.data.frame(t(melanoma_PD1_expr1)) #prepare for WGCNA  

#test genes and samples---------------------------------------------------------------------------------
gsg = goodSamplesGenes(expr, verbose = 3)
# Flagging genes and samples with too many missing values...
# ..step 1
 gsg$allOK
# [1] TRUE   
#no need to delete genes

sampleTree = hclust(dist(expr), method = "average");
jpeg('sample_cluster.jpeg',width = 1000, height = 500, units = "px", pointsize = 6,quality = 100,bg = "#e5ecff",res=100)
plot(sampleTree,sub="", xlab="")
dev.off()
#have been removed batch effect so no need to remove samples

#Choose best power-------------------------------------------------------------------------------------------------------------
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(expr, powerVector = powers, verbose = 5)

# Plot
jpeg('soft_threshold_power.jpeg',width = 1000, height = 500, units = "px", pointsize = 10,quality = 100,bg = "#e5ecff",res=100)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#power=12


#Generating adjacency and TOM similarity matrices based on the selected softpower--------------------------------------------------------
softPower = 12
adjacency = adjacency(expr, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM  # Turn adjacency into topological overlap
geneTree = hclust(as.dist(dissTOM), method = "average") ## Call the hierarchical clustering function
jpeg('TOM_based_gene_cluster.jpeg',width = 1500, height = 800, units = "px", pointsize = 12,quality = 100,bg = "#e5ecff",res=100)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04)
dev.off()
