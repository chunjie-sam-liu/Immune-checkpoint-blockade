#WGCNA for melanoma PD1
library(WGCNA)
library(magrittr)

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
tibble::rownames_to_column(all_expression) %>%
  dplyr::select(rowname,melanoma_PD1$Run)->melanoma_PD1_expr0

rownames(melanoma_PD1_expr0)=melanoma_PD1_expr0$rowname
melanoma_PD1_expr0=melanoma_PD1_expr0[,-1]
keep <- rowSums(melanoma_PD1_expr0>0) >= 2
melanoma_PD1_expr0 <- melanoma_PD1_expr0[keep,]

melanoma_PD1_expr = as.data.frame(t(melanoma_PD1_expr0))   

#test genes and samples---------------------------------------------------------------------------------
gsg = goodSamplesGenes(melanoma_PD1_expr, verbose = 3)
# Flagging genes and samples with too many missing values...
# ..step 1
 gsg$allOK
# [1] TRUE   
#no need to delete genes

sampleTree = hclust(dist(datExpr0), method = "average");
jpeg('sample_cluster.jpeg',width = 1000, height = 500, units = "px", pointsize = 6,quality = 100,bg = "#e5ecff",res=100)
plot(sampleTree,sub="", xlab="")
dev.off()
#have been removed batch effect so no need to remove samples

#Choose best power-------------------------------------------------------------------------------------------
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(melanoma_PD1_expr, powerVector = powers, verbose = 5)
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

softPower = 12
adjacency = adjacency(melanoma_PD1_expr, power = softPower)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM  # Turn adjacency into topological overlap
geneTree = hclust(as.dist(dissTOM), method = "average") ## Call the hierarchical clustering function
jpeg('TOM_based_gene_cluster.jpeg',width = 1500, height = 800, units = "px", pointsize = 12,quality = 100,bg = "#e5ecff",res=100)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04)
dev.off()
