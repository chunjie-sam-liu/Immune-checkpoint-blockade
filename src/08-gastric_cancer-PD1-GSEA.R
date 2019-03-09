library(dplyr)
library(magrittr)
library(fgsea)
library(reactome.db)


#fgsea in reactome
read.table("/data/liull/immune-checkpoint-blockade/different_expression/gastric_cancer/gastric_cancer_PD1_DEG.txt",header = T,as.is = TRUE) -> pretreatment_gastric_cancer_PD1_DEG
pretreatment_gastric_cancer_PD1_DEG %>%
  dplyr::select(ensembl_ID,t_statistic)->t_statistic
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,t_statistic,by.x="EnsemblId",by.y="ensembl_ID")%>%
  dplyr::select(GeneID,t_statistic)->gene_t_statistic

gene_list=gene_t_statistic$t_statistic
names(gene_list)=gene_t_statistic$GeneID

my_pathways <- reactomePathways(names(gene_list))
fgsea_reactome <- fgsea(pathways = my_pathways, 
                        stats = gene_list,
                        minSize=15,
                        maxSize=500,
                        nperm=100000)
head(fgsea_reactome[order(pval), ])
sum(fgsea_reactome[, padj < 0.05])


#plot the single significantly enriched pathway
plotEnrichment(my_pathways[['Extracellular matrix organization']],
               gene_list) + labs(title='Extracellular matrix organization')->ECM_organization
plotEnrichment(my_pathways[['Degradation of the extracellular matrix']],
               gene_list) + labs(title='Degradation of the extracellular matrix')->ECM_degradation
plotEnrichment(my_pathways[['Collagen formation']],
               gene_list) + labs(title='Collagen formation')->Collagen_formation
plotEnrichment(my_pathways[['Collagen degradation']],
               gene_list) + labs(title='Collagen degradation')->Collagen_degradation



ggsave(
  filename = 'Down_ECM_organization.pdf',
  plot = ECM_organization,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/GSEA/gastric_cancer_PD1/',
  width = 12,
  height = 8
)

plotEnrichment(my_pathways[['Vpu mediated degradation of CD4']],
               gene_list) + labs(title='Vpu mediated degradation of CD4') ->CD3_degradation
plotEnrichment(my_pathways[['Cross-presentation of soluble exogenous antigens (endosomes)']],
               gene_list) + labs(title='Cross-presentation of soluble exogenous antigens (endosomes)') -> endosomes
plotEnrichment(my_pathways[['Ubiquitin-dependent degradation of Cyclin D1']],
               gene_list) + labs(title='Ubiquitin-dependent degradation of Cyclin D1') -> Cyclin_D
plotEnrichment(my_pathways[['Regulation of Apoptosis']],
               gene_list) + labs(title='Regulation of Apoptosis') -> Apoptosis

ggsave(
  filename = 'Up_CD3_degradation.pdf',
  plot = CD3_degradation,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/GSEA/gastric_cancer_PD1/',
  width = 12,
  height = 8
)

#Plot the top 10 pathways enriched at the top and bottom of the ranked list, respectively
topPathwaysUp <- fgsea_reactome[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgsea_reactome[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(my_pathways[topPathways], gene_list, fgsea_reactome, 
              gseaParam = 0.5)->p

ggsave(
  filename = 'GSEA_gastric_cancer_top20.pdf',
  plot = p,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/GSEA/gastric_cancer_PD1/',
  width = 22,
  height = 8
)

