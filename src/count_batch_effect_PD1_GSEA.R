library(dplyr)
library(magrittr)
library(fgsea)
library(reactome.db)

#fgsea in reactome
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_all_DEG.txt",header = T,as.is = TRUE) -> melanoma_PD1_DEG
tibble::rownames_to_column(melanoma_PD1_DEG) %>%
  dplyr::select(rowname,logFC)->EnsemblID_logFC
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,EnsemblID_logFC,by.x="Ensembl_ID",by.y="rowname")%>%
  dplyr::select(GeneID,logFC)->GeneID_logFC

#order all of the gene by its logFC
gene_list=GeneID_logFC$logFC
names(gene_list)=GeneID_logFC$GeneID
ordered_gene_list <- gene_list[order(gene_list)]

my_pathways <- reactomePathways(names(ordered_gene_list))
fgsea_reactome <- fgsea(pathways = my_pathways, 
                        stats = ordered_gene_list,
                        minSize=15,
                        maxSize=500,
                        nperm=100000)
head(fgsea_reactome[order(pval), ])
sum(fgsea_reactome[, padj < 0.05])

dplyr::filter(fgsea_reactome,padj< 0.05)->sig_Reactome
data.frame(lapply(sig_Reactome,as.character), stringsAsFactors=FALSE)->sig_Reactome
write.table(sig_Reactome,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/GSEA/PD1_sig_Reactome.txt",quote = FALSE,sep="\t",row.names = FALSE,col.names = TRUE)#


topPathwaysUp <- fgsea_reactome[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgsea_reactome[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(my_pathways[topPathways], ordered_gene_list, fgsea_reactome, 
              gseaParam = 0.5)->p

ggsave(
  filename = 'PD1_GSEA_top20.pdf',
  plot = p,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/GSEA/',
  width = 15,
  height = 8
)
