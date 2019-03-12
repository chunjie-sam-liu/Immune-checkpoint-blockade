library(dplyr)
library(magrittr)
library(fgsea)
library(reactome.db)


#fgsea in reactome
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/PD1/pretreatment/melanoma_PD1_DEG.txt",header = T,as.is = TRUE) -> pretreatment_melanoma_PD1_DEG
pretreatment_melanoma_PD1_DEG %>%
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

dplyr::filter(fgsea_reactome,padj< 0.05)->sig_Reactome
data.frame(lapply(sig_Reactome,as.character), stringsAsFactors=FALSE)->sig_Reactome
write.table(sig_Reactome,"/data/liull/immune-checkpoint-blockade/GSEA/melanoma_PD1/pretreatment/sig_Reactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#

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
  filename = 'Down_Collagen_degradation.pdf',
  plot = Collagen_degradation,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/GSEA/melanoma_PD1/pretreatment',
  width = 12,
  height = 8
)



plotEnrichment(my_pathways[['Phosphorylation of CD3 and TCR zeta chains']],
               gene_list) + labs(title='Phosphorylation of CD3 and TCR zeta chains') ->Pho_CD3
plotEnrichment(my_pathways[['PD-1 signaling']],
               gene_list) + labs(title='PD-1 signaling') -> PD1_signaling
plotEnrichment(my_pathways[['Interferon gamma signaling']],
               gene_list) + labs(title='Interferon gamma signaling') -> IFNr_signaling
plotEnrichment(my_pathways[['TCR signaling']],
               gene_list) + labs(title='TCR signaling') -> TCR_signaling

ggsave(
  filename = 'Up_Pho_CD3.pdf',
  plot = Pho_CD3,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/GSEA/melanoma_PD1/pretreatment',
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
  filename = 'GSEA_melanoma_top20.pdf',
  plot = p,
  device = 'pdf',
  path = '//data/liull/immune-checkpoint-blockade/GSEA/melanoma_PD1/pretreatment',
  width = 15,
  height = 8
)
