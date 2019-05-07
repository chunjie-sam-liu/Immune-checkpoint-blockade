#RNA-splicing melanoma PD1
library (VennDiagram)

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP070710/A3SS.MATS.JCEC.txt",header = T) -> SRP070710_A3SS_JCEC
SRP070710_A3SS_JCEC %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP070710_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP094781/A3SS.MATS.JCEC.txt",header = T) -> SRP094781_A3SS_JCEC
SRP094781_A3SS_JCEC %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP094781_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP150548/A3SS.MATS.JCEC.txt",header = T) -> SRP150548_A3SS_JCEC
SRP150548_A3SS_JCEC %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP150548_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/ERP107734/A3SS.MATS.JCEC.txt",header = T) -> ERP107734_A3SS_JCEC
ERP107734_A3SS_JCEC %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->ERP107734_sig_rna

# intersect(SRP070710_sig_rna$geneSymbol,SRP094781_sig_rna$geneSymbol)%>%
#   intersect(SRP150548_sig_rna$geneSymbol)
venn.diagram(list(SRP070710=SRP070710_sig_rna$geneSymbol, SRP094781=SRP094781_sig_rna$geneSymbol,
                  SRP150548=SRP150548_sig_rna$geneSymbol),filename=NULL,fill = c("cornflowerblue", "green", "yellow"))->venn_melanoma_PD1
grid.draw(venn_melanoma_PD1)

# intersect(SRP070710_sig_rna$geneSymbol,SRP094781_sig_rna$geneSymbol)->a
# intersect(SRP070710_sig_rna$geneSymbol,SRP150548_sig_rna$geneSymbol)->b
# intersect(SRP094781_sig_rna$geneSymbol,SRP150548_sig_rna$geneSymbol)->c
# union(a,b)%>%
#   union(c)->all_genes

union(SRP070710_sig_rna$geneSymbol,SRP094781_sig_rna$geneSymbol)%>%
  union(SRP150548_sig_rna$geneSymbol)%>%
  union(ERP107734_sig_rna$geneSymbol)->all_genes

read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
data.frame(all_genes)%>%
  merge(relationship,by.x="all_genes",by.y="Symbol")->all_genes_frame

enrichPathway(gene=all_genes_frame$GeneID,pvalueCutoff=1, readable=T,organism = "human")->eReactome
dotplot(eReactome, showCategory=20)





#CTLA4
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP011540/SE.MATS.JCEC.txt",header = T) -> SRP011540_A3SS_JCEC
SRP011540_A3SS_JCEC %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP011540_sig_rna
