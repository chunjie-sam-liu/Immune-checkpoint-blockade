read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_up_ENSG.txt",header = T)->melanoma_PD1_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/PD1_down_ENSG.txt",header = T)->melanoma_PD1_down

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/up_ENSG.txt",header = T)->gastric_cancer_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/gastric_cancer/down_ENSG.txt",header = T)->gastric_cancer_down

read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/CTLA4_up_ENSG.txt",
           header = T)->melanoma_CTLA4_up
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_CTLA4/Second_Response_standard/CTLA4_down_ENSG.txt",
           header = T)->melanoma_CTLA4_down

intersect(melanoma_PD1_up$Symbol,gastric_cancer_up$Symbol)
#"LGALS17A" "UBD"      "LRG1"     "LCK"      "JAKMIP1"  "TRAT1"    "IDO1"     "HLA-DOA"  "PDCD1"    "HCP5"     "HLA-DQA1" "NELL2"
intersect(melanoma_PD1_down$Symbol,gastric_cancer_down$Symbol)
#"HPSE2"   "HIF3A"   "NSG1"    "RSPO4"   "GFRA3"   "CCDC8"   "RGS22"   "NTF3"    "ESYT3"   "MMP2"    "ROS1" "CORO2B" 
#"MS4A8"   "VEGFD"   "MAGED4"  "MAGED4B" "RCOR2"   "TM6SF2"  "ADAM33"  "SCUBE3"  "ANO2"    "MUC13"   "OLFM1"   "VCAN"    "ZNF662" 
intersect(melanoma_PD1_up$Symbol,gastric_cancer_up$Symbol) %>% intersect(melanoma_CTLA4_up$Symbol)
#"UBD"      "LCK"      "JAKMIP1"  "TRAT1"    "IDO1"     "HLA-DOA"  "PDCD1"    "HLA-DQA1"
intersect(melanoma_PD1_down$Symbol,gastric_cancer_down$Symbol) %>% intersect(melanoma_CTLA4_down$Symbol)
#"CORO2B"
