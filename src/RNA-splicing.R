#RNA-splicing melanoma PD1
library (VennDiagram)
library(magrittr)
library (VennDiagram)
#"MXE.MATS.JCEC.txt":6~13
#"A3SS.MATS.JCEC.txt","A5SS.MATS.JCEC.txt","MXE.MATS.JCEC.txt","RI.MATS.JCEC.txt","SE.MATS.JCEC.txt":6~11

#not consider skipping/inclusion/position------------------------------------------------ 
# file_names=c("A3SS.MATS.JCEC.txt","A5SS.MATS.JCEC.txt","MXE.MATS.JCEC.txt","RI.MATS.JCEC.txt","SE.MATS.JCEC.txt")
# 
# for (i in 1:5) {
#   
#   .x=file_names[i]
#   read.table(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP070710/",.x,sep = ""),header = T) %>%
#     dplyr::filter(FDR <=0.05)%>%
#     dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP070710_sig_rna
#   #print(paste("SRP070710",nrow(SRP070710_sig_rna),sep = ":"))
#   read.table(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP094781/",.x,sep = ""),header = T)  %>%
#     dplyr::filter(FDR <=0.05)%>%
#     dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP094781_sig_rna
#   #print(paste("SRP094781",nrow(SRP094781_sig_rna),sep = ":"))
#   read.table(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP150548/",.x,sep = ""),header = T)  %>%
#     dplyr::filter(FDR <=0.05)%>%
#     dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP150548_sig_rna
#   #print(paste("SRP150548",nrow(SRP150548_sig_rna),sep = ":"))
#   
#   unlist(strsplit(.x,"[.]"))[1] -> types
#   
#   venn.diagram(list(SRP070710=SRP070710_sig_rna$geneSymbol, 
#                     SRP094781=SRP094781_sig_rna$geneSymbol,
#                     SRP150548=SRP150548_sig_rna$geneSymbol),
#                filename=NULL,main=types,
#                fill = c("#2096BA", "#C5919D", "#DF6E21"),cat.dist=-0.1)->venn1
#   
#   ggsave(
#     filename = paste(types,"interaction.png",sep = "_"),plot = venn1,device = 'png',
#     path = '/data/liull/immune-checkpoint-blockade/RNA_splicing/',
#     width = 6,height = 6,units = "in",dpi = 100
#   )
#   
# }


#consider position and >0 / <0------------------------------------------------
#>= 0.05  Skipping
#<= -0.05 Inclusion
file_names=c("A3SS.MATS.JCEC.txt","A5SS.MATS.JCEC.txt","MXE.MATS.JCEC.txt","RI.MATS.JCEC.txt","SE.MATS.JCEC.txt")

for (i in 1:5) {
  
  .x=file_names[i]
  read.table(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP070710/",.x,sep = ""),header = T) %>%
    dplyr::filter(PValue <=0.05)%>%
    dplyr::filter(IncLevelDifference >= 0.05) ->SRP070710_sig_rna
  SRP070710_Splicing_Position=paste(SRP070710_sig_rna$geneSymbol,SRP070710_sig_rna[,6],
                                    SRP070710_sig_rna[,7],SRP070710_sig_rna[,8],SRP070710_sig_rna[,9],
                                    SRP070710_sig_rna[,10],SRP070710_sig_rna[,11],sep="_")
  #print(paste("SRP070710",nrow(SRP070710_sig_rna),sep = ":"))
  
  read.table(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP094781/",.x,sep = ""),header = T)  %>%
    dplyr::filter(PValue <=0.05)%>%
    dplyr::filter(IncLevelDifference >= 0.05) ->SRP094781_sig_rna
  SRP094781_Splicing_Position=paste(SRP094781_sig_rna$geneSymbol,SRP094781_sig_rna[,6],
                                    SRP094781_sig_rna[,7],SRP094781_sig_rna[,8],SRP094781_sig_rna[,9],
                                    SRP094781_sig_rna[,10],SRP094781_sig_rna[,11],sep="_")
  #print(paste("SRP094781",nrow(SRP094781_sig_rna),sep = ":"))
  
  read.table(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP150548/",.x,sep = ""),header = T)  %>%
    dplyr::filter(PValue <=0.05)%>%
    dplyr::filter(IncLevelDifference >= 0.05) ->SRP150548_sig_rna
  SRP150548_Splicing_Position=paste(SRP150548_sig_rna$geneSymbol,SRP150548_sig_rna[,6],
                                    SRP150548_sig_rna[,7],SRP150548_sig_rna[,8],SRP150548_sig_rna[,9],
                                    SRP150548_sig_rna[,10],SRP150548_sig_rna[,11],sep="_")
  #print(paste("SRP150548",nrow(SRP150548_sig_rna),sep = ":"))
  
  unlist(strsplit(.x,"[.]"))[1] -> types
  
  venn.diagram(list(SRP070710=SRP070710_Splicing_Position, 
                    SRP094781=SRP094781_Splicing_Position,
                    SRP150548=SRP150548_Splicing_Position),
               filename=NULL,main=types,
               fill = c("#2096BA", "#C5919D", "#DF6E21"),cat.dist=-0.05)->venn1
  
  ggsave(
    filename = paste("Skipping",types,"interaction.png",sep = "_"),plot = venn1,device = 'png',
    path = '/data/liull/immune-checkpoint-blockade/RNA_splicing/melanoma_PD1/',
    width = 6,height = 6,units = "in",dpi = 100
  )
  
}

# no consider position, interaction in two project ,enrichment----------------------------

#>= 0.05  Skipping
#<= -0.05 Inclusion

file_names=c("A3SS.MATS.JCEC.txt","A5SS.MATS.JCEC.txt","MXE.MATS.JCEC.txt","RI.MATS.JCEC.txt","SE.MATS.JCEC.txt")


two_genes=c()
for (i in 1:5) {
  
  .x=file_names[i]
  read.table(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP070710/",.x,sep = ""),header = T) %>%
    dplyr::filter(PValue <=0.05)%>%
    dplyr::filter(IncLevelDifference >= 0.05) ->SRP070710_sig_rna
  
  #print(paste("SRP070710",nrow(SRP070710_sig_rna),sep = ":"))
  
  read.table(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP094781/",.x,sep = ""),header = T)  %>%
    dplyr::filter(PValue <=0.05)%>%
    dplyr::filter(IncLevelDifference >= 0.05) ->SRP094781_sig_rna
  
  #print(paste("SRP094781",nrow(SRP094781_sig_rna),sep = ":"))
  
  read.table(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP150548/",.x,sep = ""),header = T)  %>%
    dplyr::filter(PValue <=0.05)%>%
    dplyr::filter(IncLevelDifference >= 0.05) ->SRP150548_sig_rna
  
  #print(paste("SRP150548",nrow(SRP150548_sig_rna),sep = ":"))
  
  unlist(strsplit(.x,"[.]"))[1] -> types
  intersect(SRP070710_sig_rna$geneSymbol,SRP094781_sig_rna$geneSymbol)->a
  intersect(SRP070710_sig_rna$geneSymbol,SRP150548_sig_rna$geneSymbol)->b
  intersect(SRP094781_sig_rna$geneSymbol,SRP150548_sig_rna$geneSymbol)->c
  union(a,b)%>%union(c)%>%union(two_genes)->two_genes
  
}


gsub("FAM173B","ATPSCKMT",two_genes)->two_genes
gsub("FBXO18","FBH1",two_genes)->two_genes
gsub("TROVE2","RO60",two_genes)->two_genes

gsub("C20orf24","RAB5IF",two_genes)->two_genes
gsub("METTL13","EEF1AKNMT",two_genes)->two_genes
gsub("RTFDC1","RTF2",two_genes)->two_genes
gsub("C7orf49","CYREN",two_genes)->two_genes


data.frame(Symbol=two_genes)->genes
readxl::read_xlsx("/data/liull/reference/All_EntrezID_Symbl_NCBI.xlsx") ->relationship

merge(genes,relationship)->genes


enrichGO(gene = genes$GeneID,OrgDb = org.Hs.eg.db,ont = "ALL",
         pAdjustMethod = "fdr",pvalueCutoff = 0.05,readable = TRUE)%>%
  as.data.frame(stringsAsFactors=FALSE)->enrich_GO
write.table(enrich_GO,"/data/liull/immune-checkpoint-blockade/RNA_splicing/melanoma_PD1/Skipping_two_project_enrichGO.txt",
            row.names = F,col.names = F,sep="\t")
#DOSE::dotplot(enrich_GO, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")


dplyr::mutate(enrich_GO,GO_term=paste(enrich_GO$ONTOLOGY,enrich_GO$Description,sep = ": "))-> enrich_GO

ggplot(enrich_GO,aes(x=GO_term,y=-log10(p.adjust))) +
  geom_bar(fill="#0C8DC4",stat="identity") +
  coord_flip()->enrich_GO_plot
ggsave(
  filename = 'Skipping_two_project_enrichGO.pdf',
  plot = enrich_GO_plot,
  device = 'pdf',
  path = "/data/liull/immune-checkpoint-blockade/RNA_splicing/melanoma_PD1/",
  width = 10,
  height = 6
)

#enrichKEGG(gene=genes$GeneID,organism="human",pvalueCutoff=0.05,pAdjustMethod = "BH")->enrich_KEGG#0



#number plot--------------------------------------------
#>= 0.05  Skipping
#<= -0.05 Inclusion

readxl::read_xlsx("/data/liull/immune-checkpoint-blockade/RNA_splicing/splicing.xlsx",skip = 1,sheet = "Skipping")%>%
  as.data.frame(stringsAsFactors=FALSE)->splicing_events
rownames(splicing_events)=splicing_events$Project
splicing_events[,-1] %>%
  t()%>%
  as.data.frame(stringsAsFactors=FALSE) ->splicing_events

splicing_plots=data.frame()
for (i in 1:5) {
  temp=data.frame(Project=colnames(splicing_events)[i],
                  Class=rownames(splicing_events),
                  Number=splicing_events[,i])
  
  splicing_plots=rbind(splicing_plots,temp)                
}


pdf(file = "/data/liull/immune-checkpoint-blockade/RNA_splicing/melanoma_PD1/Skipping_number.pdf", 6, 6)
ggplot(splicing_plots, aes(Project,Number, fill=Class))+
  geom_bar(stat='identity',position='dodge')+
  labs(title="Skipping in Response")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border=element_rect(fill = NA),
        axis.text.x = element_text(size = 7),
        plot.title = element_text(hjust = 0.5))
dev.off()



#gastric PD1----------------------------------------------------

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/ERP107734/A3SS.MATS.JCEC.txt",header = T) %>%
  dplyr::filter(FDR <=0.01)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->A3SS
A3SS[order(A3SS$FDR),]->A3SS
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/ERP107734/A5SS.MATS.JCEC.txt",header = T) %>%
  dplyr::filter(FDR <=0.01)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->A5SS
A5SS[order(A5SS$FDR),]->A5SS
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/ERP107734/MXE.MATS.JCEC.txt",header = T) %>%
  dplyr::filter(FDR <=0.01)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->MXE
MXE[order(MXE$FDR),]->MXE
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/ERP107734/RI.MATS.JCEC.txt",header = T) %>%
  dplyr::filter(FDR <=0.01)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->RI
RI[order(RI$FDR),]->RI
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/ERP107734/SE.MATS.JCEC.txt",header = T) %>%
  dplyr::filter(FDR <=0.01)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SE
SE[order(SE$FDR),]->SE

all=c(as.character(A3SS$geneSymbol),as.character(A5SS$geneSymbol),as.character(MXE$geneSymbol),as.character(RI$geneSymbol),as.character(SE$geneSymbol))
all[duplicated(all)]


read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_ENSG.txt",header = T)->up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_ENSG.txt",header = T)->down

#melanoma CTLA4----------------------------------------------------

read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP011540/A3SS.MATS.JCEC.txt",header = T) %>%
  dplyr::filter(FDR <=0.01)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->A3SS
A3SS[order(A3SS$FDR),]->A3SS
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP011540/A5SS.MATS.JCEC.txt",header = T) %>%
  dplyr::filter(FDR <=0.01)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->A5SS
A5SS[order(A5SS$FDR),]->A5SS
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP011540/MXE.MATS.JCEC.txt",header = T) %>%
  dplyr::filter(FDR <=0.01)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->MXE
MXE[order(MXE$FDR),]->MXE
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP011540/RI.MATS.JCEC.txt",header = T) %>%
  dplyr::filter(FDR <=0.01)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->RI
RI[order(RI$FDR),]->RI
read.table("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP011540/SE.MATS.JCEC.txt",header = T) %>%
  dplyr::filter(FDR <=0.01)%>%
  dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SE
SE[order(SE$FDR),]->SE

all=c(as.character(A3SS$geneSymbol),as.character(A5SS$geneSymbol),as.character(MXE$geneSymbol),as.character(RI$geneSymbol),as.character(SE$geneSymbol))
all[duplicated(all)]

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_up_ENSG.txt",header = T)->up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/CTLA4_down_ENSG.txt",header = T)->down
