library("survminer")

#

read.table("/data/liull/immune-checkpoint-blockade/different_expression/melanoma/melanoma_PD1_removed_batch_expression.txt",header = T,as.is = TRUE)->expr  #PD1 or CTLA4
cbind(rownames(expr),expr)->expr
colnames(expr)[1]="ensembl_ID"
read.table("/data/liull/reference/EntrezID_Symbl_EnsemblID_NCBI.txt",sep="\t",header = T,as.is = TRUE) ->relationship
merge(relationship,expr,by.x="EnsemblId",by.y="ensembl_ID",all=TRUE)%>%
  dplyr::filter(EnsemblId %in% expr$ensembl_ID) %>%
  dplyr::select(-EnsemblId,-Symbol)->expr2

Collagen_deg=c("1299", "1306", "1289", "1303", "1278", "1513", "1293", "1290", "1292", "1277", "1281", "4313", "1291", "4320", "1307")
ECM_deg=c("9510", "1299", "1306", "1634", "3915", "824", "147968", "1289", "1303", "1278", "1513", "4811", "1293", "2006", "1290", "1292", "9507", "5654", "1277", "1281", "4313", "1291", "4320", "649", "4326", "1307")
ECM_org=c("55214", "80781", "375790", "4323", "8751", "3371", "7450", "5175", "5802", "8076", "4236", "4016", "9510", "5034", "1299", "4054", "3910", "6678", "1306", "10516", "4239", "1634", "8985", "7837", "3915", "824", "147968", "1289", "1303", "85301", "1278", "3791", "1513", "2192", "4811", "1293", "2006", "1290", "3678", "1292", "3655", "9507", "5654", "871", "64175", "1277", "7043", "1281", "4313", "10536", "4060", "3679", "1291", "10491", "9509", "8974", "1462", "4320", "5339", "649", "5155", "4326", "1307", "4237")
Collagen_form=c("5034", "1299", "1306", "8985", "7837", "1289", "1303", "85301", "1278", "1293", "1290", "1292", "3655", "871", "64175", "1277", "1281", "10536", "1291", "10491", "9509", "8974", "5339", "649", "1307")
intersect(Collagen_deg,ECM_deg) %>% intersect(ECM_org) %>% intersect(Collagen_form)->Intersection_genes


readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Survival_time != "NA")%>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%###
  dplyr::select(Run,Survival_time,Survival_status) ->metadata

#filter gene and sample  c(Collagen_deg,ECM_deg,ECM_org,Collagen_form) or c("CD74","GZMA","CCL5")
expr2 %>%
  dplyr::filter(GeneID %in% c("3001","972","6352"))%>%
  dplyr::select(GeneID,metadata$Run)->expr3
apply(expr3[,2:ncol(expr3)],2,function(x) median(x))->mean_value
cbind(colnames(expr3),c("Mean",mean_value),rep("class",ncol(expr3)))[-1,]->temp
colnames(temp)=c("Run","expr_value","Class")
Combined_data=merge(metadata,as.data.frame(temp))

Combined_data$expr_value=as.numeric(as.character(Combined_data$expr_value))
Combined_data$Class=as.character(Combined_data$Class)

cutoff=median(Combined_data$expr_value)
for (i in 1:nrow(Combined_data)){
  if(Combined_data$expr_value[i] >= cutoff){
    print(i)
    Combined_data$Class[i]="high"
  }else{
    Combined_data$Class[i]="low"
  }
}


for (j in 1:nrow(Combined_data)) {
  if(Combined_data$Survival_status[j]=="Dead"){
    Combined_data$Survival_status[j]="2"
  }else {
    Combined_data$Survival_status[j]="1"
    }
}

Combined_data$Survival_time=as.numeric(Combined_data$Survival_time)
Combined_data$Survival_status=as.numeric(Combined_data$Survival_status)

fit <- survfit(Surv(Survival_time, Survival_status) ~ Class, data = Combined_data)
ggsurvplot(fit, data = Combined_data, pval = TRUE,risk.table = TRUE,risk.table.col = "strata")

