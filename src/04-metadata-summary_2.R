# library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(ggalluvial)


public_metadata <- readxl::read_xlsx(path = '/data/liull/immune-checkpoint-blockade/all_metadata_available_picture.xlsx', sheet = "SRA")
dbgap_metadata <- readxl::read_xlsx(path =  '/data/liull/immune-checkpoint-blockade/all_metadata_available_picture.xlsx', sheet = "dbGAP")

# RNA-seq data -------------------------------------------------------------
rbind(public_metadata,dbgap_metadata) %>%
  dplyr::filter(Library_strategy %in% c("RNA-Seq","single-cell seq")) -> RNA_seq_metadata

RNA_seq_metadata %>% 
  dplyr::group_by(SRA_Study) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(sra_str = glue::glue('{SRA_Study} ({n})')) %>% 
  dplyr::arrange(n) %>% 
  dplyr::select(-n) ->
  sra_str

RNA_seq_metadata %>% 
  dplyr::group_by(Cancer) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(ct_str = glue::glue('{`Cancer`} ({n})')) %>% 
  dplyr::arrange(n) %>% 
  dplyr::select(-n) -> ct_str

#%>% dplyr::mutate(ct_str = ifelse(`Cancer type` == 'NSCLC', ct_str, stringr::str_to_title(ct_str))) ->


RNA_seq_metadata %>% 
  dplyr::group_by(Library_strategy) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(lib_str = glue::glue('{Library_strategy} ({n})')) %>% 
  dplyr::arrange(n) %>% 
  dplyr::select(-n) ->
  lib_str

RNA_seq_metadata %>% 
  dplyr::group_by(`Anti_target`) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(ant_str = glue::glue('{Anti_target} ({n})')) %>% 
  dplyr::arrange(n) %>% 
  dplyr::select(-n) ->
  ant_str

RNA_seq_metadata %>% 
  dplyr::group_by(SRA_Study, `Cancer`, `Anti_target`, Library_strategy) %>% 
  dplyr::count() %>% 
  dplyr::left_join(sra_str, by = 'SRA_Study') %>% 
  dplyr::left_join(ct_str, by = 'Cancer') %>% 
  dplyr::left_join(lib_str, by = 'Library_strategy') %>% 
  dplyr::left_join(ant_str, by = 'Anti_target') %>% 
  ggplot(aes(y = n, axis1 = sra_str, axis2 = ct_str, axis3 = ant_str)) +
  geom_alluvium(aes(fill = lib_str), width = 1/12) +
  geom_stratum(width = 1/12, fill = "white", color = "red") +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c('sra_str' = 'Study', 'ct_str' = 'Cancer', 'ant_str' = 'Antibody'), expand = expand_scale(mult = c(0.18, 0.26))) +
  scale_y_continuous(expand = expand_scale(mult = c(0.01, 0.02))) +
  scale_fill_brewer(palette = "Set1") +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    
    panel.background = element_rect(fill = NA, colour = NA),
    
    legend.position = 'top'
  ) +
  guides(
    fill = guide_legend(
      title = 'Library Strategy'
    )
  ) ->
  plot_data_summary

ggsave(
  filename = '02-phase1-RNA-seq-data-summary.pdf',
  plot = plot_data_summary,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/',
  width = 6,
  height = 6.8
)



# WES_WGS medadata ----------------------------------------------------------
rbind(public_metadata,dbgap_metadata) %>%
  dplyr::filter(Library_strategy %in% c("WES","WGS")) -> WES_WGS_metadata


WES_WGS_metadata %>% 
  dplyr::group_by(SRA_Study) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(acc = glue::glue('{SRA_Study} ({n})')) %>% 
  dplyr::select(-n) ->
  acc

WES_WGS_metadata %>% 
  dplyr::group_by(Cancer) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(cc = glue::glue('{Cancer} ({n})')) %>% 
  dplyr::select(-n) ->
  cc

WES_WGS_metadata %>% 
  dplyr::group_by(Library_strategy) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(lib = glue::glue('{Library_strategy} ({n})')) %>% 
  dplyr::select(-n) ->
  lib

WES_WGS_metadata %>% 
  dplyr::group_by(`Anti_target`) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(ant = glue::glue('{`Anti_target`} ({n})')) %>% 
  dplyr::select(-n) ->
  ant

WES_WGS_metadata %>% 
  dplyr::group_by(SRA_Study, Cancer, `Anti_target`, Library_strategy) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(acc, by = 'SRA_Study') %>% 
  dplyr::left_join(cc, by = 'Cancer') %>% 
  dplyr::left_join(lib, by = 'Library_strategy') %>% 
  dplyr::left_join(ant, by = 'Anti_target') %>% 
  ggplot(aes(y = n, axis1 = acc, axis2 = cc, axis3 = ant)) +
  geom_alluvium(aes(fill = lib), width = 1/12) +
  geom_stratum(width = 1/12, fill = "white", color = "grey") +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c('acc' = 'Study', 'cc' = 'Cancer', 'ant' = 'Antibody'), expand = expand_scale(mult = c(0.24, 0.25))) +
  scale_y_continuous(expand = expand_scale(mult = c(0.01, 0.02))) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = 'WES_WGS phase2 data summary') +
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
  plot_data_WES_WGS_summary

ggsave(
  filename = '02-phase2-WES_WGS-data-summary.pdf',
  plot = plot_data_WES_WGS_summary,
  device = 'pdf',
  path = '/data/liull/immune-checkpoint-blockade/',
  width = 6,
  height = 6.8
)
