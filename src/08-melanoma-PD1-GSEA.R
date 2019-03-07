library(dplyr)
library(magrittr)
library(fgsea)
library(reactome.db)


#fgsea in reactome
read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_PD1_removed_batch_expression.txt",header = T,as.is = TRUE) -> melanoma_PD1

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata
#select response and non_response's sample id and project id
dplyr::filter(metadata,Response %in% c("CR","PR","PRCR","R")) -> response
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response

dplyr::select(melanoma_PD1,response$Run,non_response$Run)->ordered_melanoma_PD1
cbind(rownames(ordered_melanoma_PD1),ordered_melanoma_PD1)->ordered_melanoma_PD1
colnames(ordered_melanoma_PD1)[1]="EnsemblId"

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship

merge(relationship,ordered_melanoma_PD1,all=TRUE)%>%
  dplyr::filter(EnsemblId %in% ordered_melanoma_PD1$EnsemblId) %>%
  dplyr::filter(GeneID %in% relationship$GeneID)%>%
  dplyr::select(-c(1,3))->ordered_melanoma_PD1_2
rownames(ordered_melanoma_PD1_2)=ordered_melanoma_PD1_2$GeneID
ordered_melanoma_PD1_2=ordered_melanoma_PD1_2[,-1]

t_statistic=apply(ordered_melanoma_PD1_2,1,function(x) t.test(x[1:nrow(response)],x[(nrow(response)+1):length(ordered_melanoma_PD1_2)])$statistic)

my_pathways <- reactomePathways(names(t_statistic))
fgsea_reactome <- fgsea(pathways = my_pathways, 
                        stats = t_statistic,
                        minSize=15,
                        maxSize=500,
                        nperm=100000)
head(fgsea_reactome[order(pval), ])
sum(fgsea_reactome[, padj < 0.01])

#plot the single significantly enriched pathway
plotEnrichment(my_pathways[['Extracellular matrix organization']],
               t_statistic) + labs(title='Extracellular matrix organization')->ECM_organization
plotEnrichment(my_pathways[['Degradation of the extracellular matrix']],
               t_statistic) + labs(title='Degradation of the extracellular matrix')->ECM_degradation
ggsave(
  filename = 'Down_ECM_degradation.pdf',
  plot = ECM_degradation,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/GSEA/melanoma_PD1',
  width = 12,
  height = 8
)



plotEnrichment(my_pathways[['Phosphorylation of CD3 and TCR zeta chains']],
               t_statistic) + labs(title='Phosphorylation of CD3 and TCR zeta chains') ->Pho_CD3
plotEnrichment(my_pathways[['PD-1 signaling']],
               t_statistic) + labs(title='PD-1 signaling') -> PD1_signaling
plotEnrichment(my_pathways[['Interferon gamma signaling']],
               t_statistic) + labs(title='Interferon gamma signaling') -> IFNr_signaling
plotEnrichment(my_pathways[['TCR signaling']],
               t_statistic) + labs(title='TCR signaling') -> TCR_signaling

ggsave(
  filename = 'Up_TCR_signaling.pdf',
  plot = TCR_signaling,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/GSEA/melanoma_PD1',
  width = 12,
  height = 8
)

#Plot the top 10 pathways enriched at the top and bottom of the ranked list, respectively
topPathwaysUp <- fgsea_reactome[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgsea_reactome[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(my_pathways[topPathways], t_statistic, fgsea_reactome, 
              gseaParam = 0.5)->p

ggsave(
  filename = 'GSEA_melanoma_top20.pdf',
  plot = p,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/GSEA/melanoma_PD1',
  width = 12,
  height = 8
)
