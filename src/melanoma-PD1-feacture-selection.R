library(magrittr) # pipe
library(purrr)
library(DESeq2) # rse Ranged SummarizedExperiments
library(mlr) # machine learning wrapper
library(caret) # machine learning solution packages
library(glmnet) # elastic net and lasso general linear model
library(doParallel) # multiple parallel run elastic net
library(mRMRe) # mutual information
library(PCAtools) # PCA plot datasets
library(VennDiagram) # veendiagram
library(propagate) # big matrix correlation
library(biglasso) # big

setwd("C:/Users/000/Desktop/ICB_expr/PD1/")

readxl::read_excel("C:/Users/000/Desktop/ICB_expr/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  #dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response,Cancer) %>%
  dplyr::filter(Response !="NE")->metadata
read.table("C:/Users/000/Desktop/ICB_expr/melanoma_PD1_pretreatment_Symbol_count_expr.txt",sep="\t",header = T,as.is = TRUE) %>%
  tibble::rownames_to_column()->melanoma_PD1_count_expr
read.table("C:/Users/000/Desktop/ICB_expr/gastric_cancer_PD1_pretreatment_Symbol_count_expr.txt",sep="\t",header = T,as.is = TRUE) %>%
  tibble::rownames_to_column()->gastric_PD1_count_expr
merge(melanoma_PD1_count_expr,gastric_PD1_count_expr) %>%
  dplyr::select(rowname,metadata$Run) ->total_PD1
rownames(total_PD1)=total_PD1$rowname
total_PD1=total_PD1[,-1]

data.frame(Response=metadata$Response,row.names = metadata$Run)->Condition
Condition$Response %>%
  gsub("CR","R",.)%>%
  gsub("PR","R",.)%>%
  gsub("SD","NR",.)%>%
  gsub("PD","NR",.)->Condition$Response
total_PD1.se <- SummarizedExperiment(assays = as.matrix(total_PD1), 
                                     colData = Condition)
readr::write_rds(x = total_PD1.se, path = 'C:/Users/000/Desktop/ICB_expr/PD1/total_PD1.se.rds.gz', compress = 'gz')

fn_se2data.frame <- function(.se) { 
  .df <- as.data.frame(t(assay(.se))) 
  .target <- colData(x = .se)[, 'Response'] 
  .target <- ordered(x = .target, levels = c('R', 'NR')) 
  .df.t <- cbind(Response = .target, .df) 
  .df.t 
}

fn_save_se2matrix <- function(.se.name, .df) {
  .filename <- glue::glue('{gsub(".se", "", .se.name)}_{ncol(.df[, -1])}.mat.tsv')
  .mat <- as.matrix(.df[, -1])
  write.table(x = .mat, file = file.path('C:/Users/000/Desktop/ICB_expr/PD1/', .filename), sep = '\t')
  # if (file.exists(.filename)) file.remove(.filename)
  # file.symlink(from = file.path('C:/Users/000/Desktop/ICB_expr/feature_selection/', .filename), to = .filename)
}
fn_load_big.matrix <- function(.se.name, .df) {
  .ds <- glue::glue('{gsub(".se", "", .se.name)}_{ncol(.df[, -1])}')
  .filename <- glue::glue('{.ds}.mat.tsv')
  .desc <- glue::glue('{.ds}.desc')
  .bin <- glue::glue('{.ds}.bin')
  .X <- if (all(file.exists(.desc, .bin))) {
    attach.big.matrix(.desc)
  } else {
    suppressWarnings(expr = file.remove(.desc, .bin))
    setupX(
      filename = .filename,
      sep = '\t',
      header = TRUE,
      has.row.names = TRUE
    )
  }
}
fn_lasso_feats <- function(.se.name, .df) {
  # save matrix
  fn_save_se2matrix(.se.name = .se.name, .df = .df)
  
  .X <- fn_load_big.matrix(.se.name = .se.name, .df = .df)
  .y <- ifelse(.df$Response == 'R', TRUE, FALSE)
  .lassofit <- cv.biglasso(
    X = .X, y = .y, penalty = 'enet',
    family = 'binomial', screen = 'SSR', alpha = 0.1,
    seed = 1234, ncores = 2
  )
  .coef <- coef(.lassofit)[which(coef(.lassofit) != 0), ][-1]
  names(sort(x = abs(.coef), decreasing = T))
}

fn_rm_high_cor <- function(.d) {
  .cor <- bigcor(x = .d, size = 2000, fun = "cor", verbose = FALSE)
  .cor <- .cor[1:nrow(.cor), 1:ncol(.cor)]
  .high_cor <- findCorrelation(x = .cor, cutoff = 0.9)
  rm(.cor)
  gc()
  if (length(.high_cor) == 0) colnames(.d) else colnames(.d[, -.high_cor])
}

fn_select_features <- function(.se.name) {
  .se <- get(x = .se.name, envir = .GlobalEnv)
  .df <- fn_se2data.frame(.se = .se)
  # 1. sd filter
  # .feats.sd <- fn_sd_median(.df = .df)
  .feats.sd <- colnames(.df)[-1]
  # 2. lasso remove zero coef
  .feats.lasso <- fn_lasso_feats(.se.name = .se.name, .df = .df[, c('Response',.feats.sd)])
  
  # 2. High Correlation
  .feats.lasso.rmcor <- fn_rm_high_cor(.d = .df[, .feats.lasso])
  
  .feats.lasso.rmcor
}

total_PD1.feats <- fn_select_features(.se.name = "total_PD1.se")
melanoma_PD1.feats <- read.table("C:/Users/000/Desktop/ICB_expr/melanoma_PD1/features.txt",header = F) %>% as.matrix() %>% as.character()
gastric_PD1.feats <- read.table("C:/Users/000/Desktop/ICB_expr/gastric/features.txt") %>% as.matrix() %>% as.character()


write.table(total_PD1.feats,"C:/Users/000/Desktop/ICB_expr/PD1/total_PD1_features.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)

fn_resample <- function(.se, .id) {
  .d <- cbind(t(assay(.se)), as.data.frame(colData(.se)[, 'Response', drop = FALSE]))
  colnames(.d) %>% gsub("-","_",.)->colnames(.d)
  .task <- makeClassifTask(
    id = .id, data = .d,
    target = 'Response', positive = 'R'
  )
  ksvm.lrn <- makeLearner(cl = 'classif.ksvm', predict.type = 'prob')
  resamdesc.inner <- makeResampleDesc(method = 'CV', iters = 10, stratify = TRUE, predict = 'both')
  resample(
    learner = ksvm.lrn,
    task = .task,
    resampling = resamdesc.inner,
    measures = mlr::auc,
    models = TRUE,
    keep.pred = TRUE
  )
}
fn_resample(.se = total_PD1.se[total_PD1.feats, ], .id = 'total_PD1')

# Resampling: cross-validation
# Measures:             auc.test   
# [Resample] iter 1:    0.9750000  
# [Resample] iter 2:    1.0000000  
# [Resample] iter 3:    0.7222222  
# [Resample] iter 4:    0.8055556  
# [Resample] iter 5:    0.7750000  
# [Resample] iter 6:    0.9722222  
# [Resample] iter 7:    0.9722222  
# [Resample] iter 8:    1.0000000  
# [Resample] iter 9:    0.8055556  
# [Resample] iter 10:   0.6666667  
# 
# 
# Aggregated Result: auc.test.mean=0.8694444
# 
# 
# Resample Result
# Task: total_PD1
# Learner: classif.ksvm
# Aggr perf: auc.test.mean=0.8694444
# Runtime: 0.55665

fn_caret_train_model <- function(.se, method = 'svmRadial') {
  .d <- fn_se2data.frame(.se = .se)
  .model <- caret::train(
    Response ~ .,
    data = .d,
    method = "knn",
    trControl = caret::trainControl("repeatedcv", number = 10, classProbs = T),
    preProcess = c("center", "scale", "YeoJohnson"),
    tuneLength = 10
  )
  .model
}



fn_predict <- function(.model, .se) {
  .d <- fn_se2data.frame(.se = .se)
  predict(.model, .d, type = "prob")
}
fn_getROC_AUC <- function(probs, true_Y) {
  probsSort <- sort(probs, decreasing = TRUE, index.return = TRUE)
  val <- unlist(probsSort$x)
  idx <- unlist(probsSort$ix)
  
  roc_y <- true_Y[idx]
  stack_x <- cumsum(roc_y == 2) / sum(roc_y == 2)
  stack_y <- cumsum(roc_y == 1) / sum(roc_y == 1)
  
  auc <- sum((stack_x[2:length(roc_y)] - stack_x[1:length(roc_y) - 1]) * stack_y[2:length(roc_y)])
  return(list(stack_x = stack_x, stack_y = stack_y, auc = auc))
}
fn_alist <- function(.c, .pred) {
  .t <- gsub("NR",2,.c) %>% gsub("R",1,.) %>% as.numeric()
  aList <- fn_getROC_AUC(probs = .pred$R, true_Y = .t)
}
fn_plot_auc <- function(.alist, .t) {
  tibble::tibble(
    x = .alist$stack_x,
    y = .alist$stack_y
  ) %>%
    ggplot(aes(x = x, y = y)) +
    geom_path() +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), size = 0.05, linetype = 2) +
    theme_bw() +
    labs(
      x = 'False Positive Rate',
      y = 'True Positive Rate',
      title = glue::glue("{.t}, AUC = {round(.alist$auc,3)}")
    )
}
fn_ensembl_predict <- function(.model, .se, .title) {
  .d.pred <- fn_predict(.model = .model, .se = .se)
  .d.pred.alist <- fn_alist(.c = .se$Response, .pred = .d.pred)
  .d.p <- fn_plot_auc(.alist = .d.pred.alist, .t = .title)
  list(predict = .d.pred, alist = .d.pred.alist, plot = .d.p)
}

smart.model <- fn_caret_train_model(.se = melanoma_PD1.se[total_PD1.feats, ])
panel.smart.test.se.auc <- fn_ensembl_predict(.model = smart.model, .se = melanoma_PD1.se[total_PD1.feats, ], .title = 'melanoma_PD1')
panel.smart.test.se.auc$plot

smart.model <- fn_caret_train_model(.se = gastric_PD1.se[total_PD1.feats, ])
panel.smart.test.se.auc <- fn_ensembl_predict(.model = smart.model, .se = gastric_PD1.se[total_PD1.feats, ], .title = 'gastric_PD1')
panel.smart.test.se.auc$plot


# fn_pca_biplot_individual <- function(.se, .title) {
#   .se.m <- assay(.se)
#   .metadata <- colData(.se)[, 'Response', drop = FALSE]
#   .p <- pca(mat = .se.m, metadata = .metadata, removeVar = 0.2)
#   .biplot <- biplot(
#     .p,
#     colby = 'Response',
#     lab = FALSE,
#     legendPosition = 'right',
#     title = glue::glue("{.title} PCA")
#   )
#   .biplot
# }
# fn_pca_biplot_individual(.se = melanoma_PD1.se, .title = 'melanoma_PD1')
# fn_pca_biplot_individual(.se = gastric_PD1.se, .title = 'gastric_PD1')

# melanoma_PD1.se <- readr::read_rds(path = "C:/Users/000/Desktop/ICB_expr/melanoma_PD1/melanoma_PD1.se.rds.gz")
# gastric_PD1.se <- readr::read_rds(path = "C:/Users/000/Desktop/ICB_expr/gastric/gastric.se.rds.gz")
# 
# fn_resample(.se = melanoma_PD1.se[total_PD1.feats, ], .id = 'melanoma_PD1')
# fn_resample(.se = gastric_PD1.se[total_PD1.feats, ], .id = 'gastric_PD1')
