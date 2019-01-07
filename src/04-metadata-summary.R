
# library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(ggalluvial)


# public data -------------------------------------------------------------

public_metadata <- readxl::read_xlsx(path = 'data/02-metadata-verified-cj.xlsx', sheet = 1)

public_metadata %>% 
  dplyr::group_by(SRA_Study) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(sra_str = glue::glue('{SRA_Study} ({n})')) %>% 
  dplyr::arrange(n) %>% 
  dplyr::select(-n) ->
  sra_str

public_metadata %>% 
  dplyr::group_by(`Cancer type`) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(ct_str = glue::glue('{`Cancer type`} ({n})')) %>% 
  dplyr::arrange(n) %>% 
  dplyr::select(-n) %>% 
  dplyr::mutate(ct_str = ifelse(`Cancer type` == 'NSCLC', ct_str, stringr::str_to_title(ct_str))) ->
  ct_str

public_metadata %>% 
  dplyr::group_by(Library_strategy) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(lib_str = glue::glue('{Library_strategy} ({n})')) %>% 
  dplyr::arrange(n) %>% 
  dplyr::select(-n) ->
  lib_str

public_metadata %>% 
  dplyr::group_by(`anti-target`) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(ant_str = glue::glue('{`anti-target`} ({n})')) %>% 
  dplyr::arrange(n) %>% 
  dplyr::select(-n) ->
  ant_str

public_metadata %>% 
  dplyr::group_by(SRA_Study, `Cancer type`, `anti-target`, Library_strategy) %>% 
  dplyr::count() %>% 
  dplyr::left_join(sra_str, by = 'SRA_Study') %>% 
  dplyr::left_join(ct_str, by = 'Cancer type') %>% 
  dplyr::left_join(lib_str, by = 'Library_strategy') %>% 
  dplyr::left_join(ant_str, by = 'anti-target') %>% 
  ggplot(aes(y = n, axis1 = sra_str, axis2 = ct_str, axis3 = ant_str)) +
  geom_alluvium(aes(fill = lib_str), width = 1/12) +
  geom_stratum(width = 1/12, fill = "white", color = "grey") +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c('sra_str' = 'Study', 'ct_str' = 'Cancer Type', 'ant_str' = 'Antibody'), expand = expand_scale(mult = c(0.18, 0.26))) +
  scale_y_continuous(expand = expand_scale(mult = c(0.01, 0.02))) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = 'public phase1 data summary') +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    
    panel.background = element_rect(fill = NA, colour = NA),
    
    legend.position = 'top',
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(
    fill = guide_legend(
      title = 'Library Strategy'
    )
  ) ->
  plot_public_data_summary
  
ggsave(
  filename = '01-phase1-pulic-data-summary.pdf',
  plot = plot_public_data_summary,
  device = 'pdf',
  path = 'figures',
  width = 6,
  height = 6.8
)



# dbGaP medadata ----------------------------------------------------------

dbgap_metadata <- readxl::read_xlsx(path = 'data/03-dbGAP-matedata-cj.xlsx', sheet = 1)

dbgap_metadata %>% 
  dplyr::group_by(Second_accession) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(acc = glue::glue('{Second_accession} ({n})')) %>% 
  dplyr::select(-n) ->
  acc

dbgap_metadata %>% 
  dplyr::group_by(Cancer) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(cc = glue::glue('{Cancer} ({n})')) %>% 
  dplyr::select(-n) ->
  cc

dbgap_metadata %>% 
  dplyr::group_by(Library_strategy) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(lib = glue::glue('{Library_strategy} ({n})')) %>% 
  dplyr::select(-n) ->
  lib

dbgap_metadata %>% 
  dplyr::group_by(`anti-target`) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(ant = glue::glue('{`anti-target`} ({n})')) %>% 
  dplyr::select(-n) ->
  ant

dbgap_metadata %>% 
  dplyr::group_by(Second_accession, Cancer, `anti-target`, Library_strategy) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(acc, by = 'Second_accession') %>% 
  dplyr::left_join(cc, by = 'Cancer') %>% 
  dplyr::left_join(lib, by = 'Library_strategy') %>% 
  dplyr::left_join(ant, by = 'anti-target') %>% 
  ggplot(aes(y = n, axis1 = acc, axis2 = cc, axis3 = ant)) +
  geom_alluvium(aes(fill = lib), width = 1/12) +
  geom_stratum(width = 1/12, fill = "white", color = "grey") +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c('acc' = 'dbGaP', 'cc' = 'Cancer Type', 'ant' = 'Antibody'), expand = expand_scale(mult = c(0.24, 0.25))) +
  scale_y_continuous(expand = expand_scale(mult = c(0.01, 0.02))) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = 'dbGaP phase2 data summary') +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    
    panel.background = element_rect(fill = NA, colour = NA),
    
    legend.position = 'top',
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(
    fill = guide_legend(
      title = 'Library Strategy'
    )
  ) ->
  plot_data_dbgap_summary

ggsave(
  filename = '02-phase2-dbgap-data-summary.pdf',
  plot = plot_data_dbgap_summary,
  device = 'pdf',
  path = 'figures',
  width = 6,
  height = 6.8
)
  