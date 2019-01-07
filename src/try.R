path_data="/data/liull/data/FPKM/"
expression_files2 <- list.files(path = path_data, pattern = 'NA_number.txt', recursive = TRUE,full.names = TRUE)
expression_files2[1] %>% 
  read.table(sep="\t",header = T,as.is = TRUE) -> data
data$non_NA_sample_number_per_gene=as.numeric(data$non_NA_sample_number_per_gene)

data=data[1:(dim(data)[1]-1),1:(dim(data)[2]-1)]#delete NA number col and row
row.names(data)=data[,1]
data=data[,-1]#make rownames
dim(data)      
data[is.na(data)] <- 0

del=c(50,51,52,53,54,57,58,59,60,61,62,65,66,67)
data2=data[,-del]#delete the wrong sample with too little gene mapped

Sum=apply(data2,1,function(x) sum(x,na.rm = TRUE))
IDs=which(Sum==0)
data3=data2[-IDs,]#delete the gene has all 0.000 depression


#batch1=rep(1,47)
#batch2=rep(2,23)
#batch3=rep(3,28)
#batch=c(batch1,batch2,batch3)
#data3=as.matrix(data3)
#combat_edata = ComBat(dat=data3, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

readxl::read_excel("/data/liull/data/FPKM/04-all-metadata.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::select(Run,Response) ->metadata
#筛选出黑色素瘤的RNA-seq的anti-PD1的数据的应答信息
CR=dplyr::filter(metadata,Response=="CR")
PR=dplyr::filter(metadata,Response=="PR")
R=dplyr::filter(metadata,Response=="R")
rbind(CR,PR) %>%
  rbind(R) ->response
SD=dplyr::filter(metadata,Response=="SD")
PD=dplyr::filter(metadata,Response=="PD")
NR=dplyr::filter(metadata,Response=="NR")
rbind(SD,PD) %>%
  rbind(NR) ->non_response
#select response and non_response

#（1）去批次前找差异
dplyr::select(data3,one_of(non_response$Run)) ->non_response_expression
x<-c(paste("control_",1:ncol(non_response_expression),sep=""))

dplyr::select(data3,one_of(response$Run)) ->response_expression
y<-c(paste("case_",1:ncol(response_expression),sep=""))

all_expression=cbind(non_response_expression,response_expression)
myfactors <- data.frame(Tissue = c(rep("control",ncol(non_response_expression)),rep("case",ncol(response_expression))),TissueRun = c(x,y))

mydata <- readData(data = all_expression,factors = myfactors)
mynoiseq <- noiseqbio(mydata, k = 0.5, norm = "n", factor = "Tissue",lc = 0, r = 20, adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345,filter = 1,cv.cutoff = 100, cpm = 0)
noiseq_result<-mynoiseq@results[[1]]
noiseq_degene<-subset(noiseq_result,noiseq_result$prob>=0.75)
dim(noiseq_degene)
#(2)combat_edata,去批次后找差异
dplyr::select(combat_edata,one_of(non_response$Run)) ->non_response_expression
x<-c(paste("control_",1:ncol(non_response_expression),sep=""))

dplyr::select(combat_edata,one_of(response$Run)) ->response_expression
y<-c(paste("case_",1:ncol(response_expression),sep=""))

all_expression=cbind(non_response_expression,response_expression)
for (i in 1:dim(all_expression)[1]) {
  replace_IDs=which(all_expression[i,]<0.01)
  all_expression[i,replace_IDs]=0.01
}
myfactors <- data.frame(Tissue = c(rep("control",ncol(non_response_expression)),rep("case",ncol(response_expression))),TissueRun = c(x,y))

mydata <- readData(data = all_expression,factors = myfactors)
mynoiseq <- noiseqbio(mydata, k = 0.5, norm = "n", factor = "Tissue",lc = 0, r = 20, adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345,filter = 1,cv.cutoff = 100, cpm = 0)
noiseq_result<-mynoiseq@results[[1]]
noiseq_degene<-subset(noiseq_result,noiseq_result$prob>=0.75)
dim(noiseq_degene)
