library(NOISeq)
library(magrittr)
library(readr)
library(readxl)
library(dplyr)
path_data <- '/data/liull/data/FPKM'
expression_files <- list.files(path = path_data, pattern = 'expression_deleted_batch.txt', recursive = TRUE, full.names = TRUE)
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

expression_files %>%
  head(1) %>%
  read.table(sep="\t",header = T,as.is = TRUE) ->data1

dplyr::select(data1,one_of(non_response$Run)) ->non_response_expression
x<-c(paste("control_",1:ncol(non_response_expression),sep=""))

dplyr::select(data1,one_of(response$Run)) ->response_expression
y<-c(paste("case_",1:ncol(response_expression),sep=""))

all_expression=cbind(non_response_expression,response_expression)
myfactors <- data.frame(Tissue = c(rep("control",ncol(non_response_expression)),rep("case",ncol(response_expression))),TissueRun = c(x,y))


for (i in 1:dim(all_expression)[1]) {
  replace_IDs=which(all_expression[i,]<0.01)
  all_expression[i,replace_IDs]=0.01
}


mydata <- readData(data = all_expression,factors = myfactors)
mynoiseq <- noiseqbio(mydata, k = 0.5, norm = "n", factor = "Tissue",lc = 0, r = 20, adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345,filter = 1,cv.cutoff = 100, cpm = 0)
noiseq_result<-mynoiseq@results[[1]]
noiseq_degene<-subset(noiseq_result,noiseq_result$prob>=0.75)
noiseq_up<-subset(noiseq_degene,noiseq_degene$log2FC>=1)
noiseq_down<-subset(noiseq_degene,noiseq_degene$log2FC<=(-1))
