
# library -----------------------------------------------------------------

library(magrittr)

# path --------------------------------------------------------------------

path_data <- '/home/liucj/data/data/immune-checkpoint-blockade'
path_fastq <- file.path(path_data, 'fastq')

# load data ---------------------------------------------------------------

metadata <- readxl::read_xlsx(path = file.path(path_data, '02-metadata-verified.xlsx')) %>% 
  dplyr::select(lib = `Library_strategy`, run = `Run`)

path_projects <- list.dirs(path = path_fastq, recursive = F)

path_projects %>% 
  purrr::map(
    .f = function(.x) {
      # create wes and rna-seq dir
      .path_wes <- file.path(.x,'WES')
      if (!dir.exists(.path_wes)) dir.create(.path_wes)
      .path_rna_seq <- file.path(.x, 'RNA-Seq')
      if (!dir.exists(.path_rna_seq)) dir.create(.path_rna_seq)
      
      
      .fastq_gz <- list.files(path = .x, pattern = 'fastq.gz')
      .fastq <- .fastq_gz %>% stringr::str_split(pattern = "[^[:alnum:]]+", simplify = TRUE) %>% .[, 1]
      
      metadata %>% dplyr::filter(run %in% .fastq) -> .fastq_meta 
      
      .fastq_meta %>% dplyr::filter(lib == 'RNA-Seq') %>% dplyr::pull(run) -> .fastq_rna_seq
      if (length(.fastq_rna_seq) != 0) system(glue::glue("cd {.x} && mv {paste(.fastq_rna_seq, '*', sep = '', collapse = ' ')} {.path_rna_seq}"))
      
      .fastq_meta %>% dplyr::filter(lib == 'WES') %>% dplyr::pull(run) -> .fastq_wes
      if (length(.fastq_wes) != 0) system(glue::glue("cd {.x} && mv {paste(.fastq_wes, '*', sep = '', collapse = ' ')} {.path_wes}"))
    }
  )
