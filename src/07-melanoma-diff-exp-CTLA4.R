library(magrittr)
library(readr)
library(readxl)
library(dplyr)
library(sva)

#筛选出黑色素瘤的RNA-seq的anti-CTLA4的数据的应答信息
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-CTLA4") %>%
  dplyr::select(SRA_Study,Run,Response) ->metadata
#dim(metadata)
#[1] 6 3
