#RNA-splicing melanoma PD1
library (VennDiagram)

#"MXE.MATS.JCEC.txt":6~13

#"A3SS.MATS.JCEC.txt","A5SS.MATS.JCEC.txt","RI.MATS.JCEC.txt","SE.MATS.JCEC.txt":6~11

.x="SE.MATS.JCEC.txt"
      
      read.table(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP070710/",.x,sep = ""),header = T) %>%
        dplyr::filter(PValue <=0.05)%>%
        dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP070710_sig_rna
      SRP070710_Splicing_Position=paste(SRP070710_sig_rna$geneSymbol,SRP070710_sig_rna[,6],
                                        SRP070710_sig_rna[,7],SRP070710_sig_rna[,8],SRP070710_sig_rna[,9],
                                        SRP070710_sig_rna[,10],SRP070710_sig_rna[,11],sep="_")
      SRP070710_sig_rna %>% dplyr::mutate(SRP070710_Splicing_Position)->SRP070710_sig_rna
      
      read.table(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP094781/",.x,sep = ""),header = T)  %>%
        dplyr::filter(PValue <=0.05)%>%
        dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP094781_sig_rna
      SRP094781_Splicing_Position=paste(SRP094781_sig_rna$geneSymbol,SRP094781_sig_rna[,6],
                                        SRP094781_sig_rna[,7],SRP094781_sig_rna[,8],SRP094781_sig_rna[,9],
                                        SRP094781_sig_rna[,10],SRP094781_sig_rna[,11],sep="_")
      SRP094781_sig_rna %>% dplyr::mutate(SRP094781_Splicing_Position)->SRP094781_sig_rna
      
      
      read.table(paste("/data/liull/immune-checkpoint-blockade/RNA_splicing/SRP150548/",.x,sep = ""),header = T)  %>%
        dplyr::filter(PValue <=0.05)%>%
        dplyr::filter(abs(IncLevelDifference) >= 0.05) ->SRP150548_sig_rna
      SRP150548_Splicing_Position=paste(SRP150548_sig_rna$geneSymbol,SRP150548_sig_rna[,6],
                                        SRP150548_sig_rna[,7],SRP150548_sig_rna[,8],SRP150548_sig_rna[,9],
                                        SRP150548_sig_rna[,10],SRP150548_sig_rna[,11],sep="_")
      SRP150548_sig_rna %>% dplyr::mutate(SRP150548_Splicing_Position)->SRP150548_sig_rna
      
#3'-----------------------------
#"OSGEPL1"  "TOR1AIP1" "RHOT2"    "RPAIN"   
# [1] "OSGEPL1_189746659_189747331_189746659_189747168_189750549_189750656"
      #protein:Involved in mitochondrial genome maintenance.
# [2] "TOR1AIP1_179889309_179889369_179889312_179889369_179884691_179884769"
      #This gene encodes a type 2 integral membrane protein that binds A- and B-type lamins

#5'-----------------------------

#"ACTL6A"  "HYI"     "ARL6IP4" "ULK3" 
# [1] "ACTL6A_179562922_179563117_179562922_179562973_179569823_179569900"  
      #This gene encodes a family member of actin-related proteins (ARPs)
# [2] "HYI_43453385_43453497_43453412_43453497_43452204_43452319" 
      #may be involved in carbohydrate transport and metabolism
# [3] "ULK3_74837750_74837798_74837756_74837798_74837368_74837435"
      #protein:Also acts as a regulator of autophagy: following cellular senescence, able to induce autophagy        


#MXE-----------------------------

#[1] "MCUB"    "NDUFS1"  "PPOX"    "ZMYND8"  "SLC37A3" "GUSB"    "GIPC1"   "CERS4"   "SUN1"    "HDAC8"   "RALY"    "ATG4B"   "PPIL3"

# "MCUB_109659010_109659086_109660194_109660365_109560204_109560436_109664289_109664394" 
  #mitochondrial calcium uniporter dominant negative beta subunit
# "ATG4B_241651263_241651335_241653511_241653610_241651009_241651111_241654545_241654647"
      #autophagy related 4B cysteine peptidase 

#RI-----------------------------

#0
#0
#SE-----------------------------

#95
# [1] "SPAG9_51031680_51031722_51021157_51021365_51041500_51041651" #This protein also binds to kinesin-1 and may be involved in microtubule-based membrane transport.       
# [2] "IQCK_19763847_19764112_19733697_19733827_19788837_19788922" #The IQ motif serves as a binding site for different EF-hand proteins such as calmodulin        
# [3] "CAST_96740744_96740783_96740037_96740118_96741265_96741317" #The protein encoded by this gene is an endogenous calpain (calcium-dependent cysteine protease) inhibitor#        
# [4] "COPZ2_48033210_48033302_48029124_48029176_48033862_48033872"#COPI vesicles function in the retrograde Golgi-to-ER transport of dilysine-tagged proteins        
# [5] "POLR2J4_44014605_44014783_44013592_44013679_44016433_44016523" #RNA polymerase II subunit J4, pseudogene      
# [6] "DRAM2_111139500_111139666_111137522_111137586_111140037_111140090"#The protein encoded by this gene binds microtubule-associated protein 1 light chain 3 and is required for autophagy.   
# [7] "DNM2_10808568_10808580_10805915_10805967_10812263_10812377"#dynamin 2,They have been implicated in cell processes such as endocytosis and cell motility, and in alterations of the membrane that accompany certain activities such as bone resorption by osteoclasts         
# [8] "GEMIN7_45079906_45080029_45079257_45079279_45090106_45090537"       
# [9] "AMDHD2_2527572_2527615_2520983_2521123_2527772_2527836"             
# [10] "DNAJC19_180986942_180987022_180985409_180985996_180988022_180988096"#The protein encoded by this gene is thought to be part of a complex involved in the ATP-dependent transport of transit peptide-containing proteins from the inner cell membrane to the mitochondrial matrix. 
# [11] "MFF_227347225_227347384_227340291_227340380_227355676_227355761" #a nuclear gene encoding a protein that functions in mitochondrial and peroxisomal fission.   
# [12] "MFF_227352513_227352573_227340291_227340380_227355676_227355761"    
# [13] "EAF1_15432086_15432223_15429912_15430007_15434347_15434538"         
# [14] "PLEC_143922503_143925884_143915146_143922395_143926783_143926882"  #capable of interlinking different elements of the cytoskeleton 
# [15] "APLP2_130123611_130123779_130122304_130122513_130126699_130126830"  
# [16] "SENP6_75640683_75640704_75634706_75634811_75647730_75647801"        
# [17] "SEC31A_82861630_82861708_82857688_82857764_82862533_82862572" #is a component of the outer layer of the coat protein complex II (COPII).      
# [18] "METTL26_635280_635340_634663_634797_636093_636278"                  
# [19] "ZNF207_32366664_32366757_32365329_32365487_32367771_32368014"       
# [20] "NCOA1_24768078_24768135_24762686_24762776_24768223_24768543"        
# [21] "DMAC2_41433271_41433434_41431318_41432408_41433536_41433673"        
# [22] "WARS_100374135_100374244_100369086_100369258_100375282_100375347"   
# [23] "WARS_100374135_100374244_100369086_100369258_100376259_100376334"   
# [24] "GMIP_19634506_19634703_19633802_19634190_19634791_19634929"         
# [25] "PPP1R7_241153482_241153604_241150466_241150547_241157806_241157862" 
# [26] "KAT5_65713347_65713503_65712921_65713058_65713592_65713664" #This protein is a histone acetylase that has a role in DNA repair and apoptosis and is thought to play an important role in signal transduction. 

# venn.diagram(list(SRP070710=SRP070710_sig_rna$geneSymbol, SRP094781=SRP094781_sig_rna$geneSymbol,
#                   SRP150548=SRP150548_sig_rna$geneSymbol),filename=NULL,fill = c("cornflowerblue", "green", "yellow"))->venn_melanoma_PD1
# grid.draw(venn_melanoma_PD1)

intersect(SRP070710_sig_rna$geneSymbol,SRP094781_sig_rna$geneSymbol)%>%
  intersect(SRP150548_sig_rna$geneSymbol)

intersect(SRP070710_sig_rna$SRP070710_Splicing_Position,SRP094781_sig_rna$SRP094781_Splicing_Position)%>%
  intersect(SRP150548_sig_rna$SRP150548_Splicing_Position)->interaction_splicing_gene



#enrichment
# read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
# data.frame(interaction_splicing_gene)%>%
#   merge(relationship,by.x="interaction_splicing_gene",by.y="Symbol")->all_genes_frame



