library(dplyr)
library(magrittr)
library(fgsea)
library(reactome.db)


#fgsea in reactome
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/CTLA4/melanoma_CTLA4_DEG.txt",header = T,as.is = TRUE) -> pretreatment_melanoma_CTLA4_DEG
pretreatment_melanoma_CTLA4_DEG %>%
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


#Plot the top 10 pathways enriched at the top and bottom of the ranked list, respectively
topPathwaysUp <- fgsea_reactome[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgsea_reactome[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(my_pathways[topPathways], gene_list, fgsea_reactome, 
              gseaParam = 0.5)->p

ggsave(
  filename = 'GSEA_melanoma_CTLA4_top20.pdf',
  plot = p,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/GSEA/melanoma_CTLA4',
  width = 22,
  height = 8
)



