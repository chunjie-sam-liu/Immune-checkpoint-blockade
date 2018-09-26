# the data is from cell paper 10.1016/j.cell.2017.01.010
# GSE78220


# library -----------------------------------------------------------------

library(magrittr)
library(GEOquery)
library(rvest)

# data path ---------------------------------------------------------------

path_data <- "/data/liucj/data/immune-checkpoint-blockade"

path_gse <- file.path(path_data, "GSE78220")

if (!file.exists(path_gse)) dir.create(path_gse)
# get data ----------------------------------------------------------------

GSE78220 <- getGEO(GEO = 'GSE78220', destdir = path_gse, GSEMatrix = TRUE)

GSE78220[[1]] %>% 
  phenoData() %>% 
  pData() ->
  gse_pdata

gse_pdata %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate(
    srr = purrr::map(
      .x = relation.1,
      .f = function(.x) {
        .srx <- gsub(pattern = "SRA: ", replacement = "", x = .x)
        .srx_html <- xml2::read_html(x = .srx)
        .srx_html %>% 
          rvest::html_nodes('td a') %>% 
          rvest::html_text()
      }
    )
  ) %>% 
  tidyr::unnest() ->
  gse_pdata_srr


# save annotation to tsv --------------------------------------------------

gse_pdata_srr %>% readr::write_tsv(path = file.path(path_gse, '01-GSE78220-annotation-file.tsv'))


# download raw data -------------------------------------------------------


gse_pdata_srr$srr %>% paste0(collapse = ',') -> srrs
cmd <- glue::glue("bash /data/liucj/github/useful-scripts/downloadsra.sh {srrs} {path_gse}")
system(cmd)
