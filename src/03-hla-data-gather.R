
# library -----------------------------------------------------------------

library(magrittr)


# path --------------------------------------------------------------------

path_hla <- '/project/miaoyr/project/01-hla-typing/'

# hla tsv -----------------------------------------------------------------

hla_list <- '/project/miaoyr/project/01-hla-typing/cromwell-executions/HLA_typing/a331516f-e2ee-47af-a796-72175ffd795d/call-gatherFile/execution/Optitype.result.tsv.list'

hla_tsvs <- readr::read_lines(file = hla_list)

hla_tsvs %>%
  purrr::map(
    .f = function(.x) {
      .srr <- basename(.x) %>% stringr::str_replace(pattern = '_result.tsv', replacement = '')
      
      .x %>% 
        readr::read_tsv() %>% 
        dplyr::select(-X1) %>% 
        tibble::add_column(sample = .srr, .before = 1)
    }
  ) %>% 
  dplyr::bind_rows() ->
  hla_rnaseq_exome_phase1


# save data ---------------------------------------------------------------

hla_rnaseq_exome_phase1 %>% readr::write_tsv(path = file.path(path_hla, '01-hla-rnaseq-exome-phase1.tsv'))
