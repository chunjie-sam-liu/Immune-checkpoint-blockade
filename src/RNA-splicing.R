#RNA-splicing melanoma PD1
library (VennDiagram)


#3'-----------------------------
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP070710/A3SS.MATS.JCEC.txt",header = T) %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP070710_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP094781/A3SS.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP094781_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP150548/A3SS.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP150548_sig_rna
#"OSGEPL1"  "TOR1AIP1" "RHOT2"    "RPAIN"

#5'-----------------------------
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP070710/A5SS.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP070710_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP094781/A5SS.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP094781_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP150548/A5SS.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP150548_sig_rna
#"ACTL6A"  "HYI"     "ARL6IP4" "ULK3" 

#MXE-----------------------------
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP070710/MXE.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP070710_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP094781/MXE.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP094781_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP150548/MXE.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP150548_sig_rna
#[1] "MCUB"    "NDUFS1"  "PPOX"    "ZMYND8"  "SLC37A3" "GUSB"    "GIPC1"   "CERS4"   "SUN1"    "HDAC8"   "RALY"    "ATG4B"   "PPIL3"

#RI-----------------------------
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP070710/RI.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP070710_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP094781/RI.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP094781_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP150548/RI.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP150548_sig_rna
#0

#SE-----------------------------
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP070710/SE.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP070710_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP094781/SE.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP094781_sig_rna

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP150548/SE.MATS.JCEC.txt",header = T)  %>%
  dplyr::filter(PValue <=0.05)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP150548_sig_rna
#95

# venn.diagram(list(SRP070710=SRP070710_sig_rna$geneSymbol, SRP094781=SRP094781_sig_rna$geneSymbol,
#                   SRP150548=SRP150548_sig_rna$geneSymbol),filename=NULL,fill = c("cornflowerblue", "green", "yellow"))->venn_melanoma_PD1
# grid.draw(venn_melanoma_PD1)

union(SRP070710_sig_rna$geneSymbol,SRP094781_sig_rna$geneSymbol)%>%
  union(SRP150548_sig_rna$geneSymbol)->all_genes
intersect(SRP070710_sig_rna$geneSymbol,SRP094781_sig_rna$geneSymbol)%>%
  intersect(SRP150548_sig_rna$geneSymbol)->interaction_splicing_gene


#DEG genes
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_down_ENSG.txt",
           header = T)->melanoma_PD1_down
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_up_ENSG.txt",
           header = T)->melanoma_PD1_up

union(melanoma_PD1_down$Symbol,melanoma_PD1_up$Symbol) -> all_DEG
intersect(all_DEG,all_genes)


#enrichment
# read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
# data.frame(interaction_splicing_gene)%>%
#   merge(relationship,by.x="interaction_splicing_gene",by.y="Symbol")->all_genes_frame



