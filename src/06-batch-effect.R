source("https://bioconductor.org/biocLite.R")
biocLite("sva")
library(sva)
library(readr)
library(magrittr)
library(dplyr)
path_data="/data/liull/data/FPKM/"
expression_files <- list.files(path = path_data, pattern = 'expression.txt', recursive = TRUE,full.names = TRUE)

for (fi in 1:2) {

read.table(expression_files[fi],sep="\t",header = T,as.is = TRUE) -> nc_edata1

nc_edata1=nc_edata1[-dim(nc_edata1)[1],]

non_NA_sample_number_per_gene=NULL
for (i in 1:dim(nc_edata1)[1]){non_NA_sample_number_per_gene[i]=112-sum(is.na(nc_edata1[i,]))}
nc_edata2=cbind(nc_edata1,non_NA_sample_number_per_gene)
non_NA_gene_number_per_sample=NULL
non_NA_gene_number_per_sample[1]="non_NA_gene_number_per_sample"
for (i in 2:(dim(nc_edata1)[2])){non_NA_gene_number_per_sample[i]=dim(nc_edata1)[1]-sum(is.na(nc_edata1[,i]))}
non_NA_gene_number_per_sample[dim(nc_edata2)[2]]=sum(as.numeric(non_NA_gene_number_per_sample[2:dim(nc_edata1)[2]]))

nc_edata2=rbind(nc_edata2,non_NA_gene_number_per_sample)
name=strsplit(expression_files[fi],"[/]")[[1]][7]
write.table(nc_edata2,paste(path_data,name,"_non_NA_number.txt",sep = ""),quote = FALSE,sep="\t",row.names = FALSE)

}


library(purrr)
expression_files2 <- list.files(path = path_data, pattern = 'NA_number.txt', recursive = TRUE,full.names = TRUE)
expression_files2 %>% 
  purrr::map(
    .f = function(.x) {
  read.table(.x,sep="\t",header = T,as.is = TRUE) -> data
data$non_NA_sample_number_per_gene=as.numeric(data$non_NA_sample_number_per_gene)
dplyr::filter(data,non_NA_sample_number_per_gene>=3*(dim(data)[2]-2)/4) ->data2
#filter the gene has too many NA

data2=data2[1:(dim(data2)[1]-1),1:(dim(data2)[2]-1)]#delete NA number col and row

row.names(data2)=data2[,1]
data2=data2[,-1]#make rownames
dim(data2)

del=c(50,51,52,53,54,57,58,59,60,61,62,65,66,67)
data2=data2[,-del]#delete the wrong sample with too little gene mapped



Sum=apply(data2,1,function(x) sum(x,na.rm = TRUE))
IDs=which(Sum==0)
data3=data2[-IDs,]#delete the gene has all 0.000 depression

a=NULL
for(i in 1:dim(data3)[2]) {
  a <- mean(data3[, i], na.rm = T)
  data3[is.na(data3[, i]), i] <- a
}#take NA to mean of its sample



batch1=rep(1,47)
batch2=rep(2,23)
batch3=rep(3,28)
batch=c(batch1,batch2,batch3)
data3=as.matrix(data3)
combat_edata = ComBat(dat=data3, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

output_name=strsplit(strsplit(.x,"[/]")[[1]][7],"[_]")[[1]][1]
write.table(combat_edata,paste(path_data,output_name,"_expression_deleted_batch.txt",sep = ""),quote = FALSE,sep="\t",row.names = TRUE)
 }
)
