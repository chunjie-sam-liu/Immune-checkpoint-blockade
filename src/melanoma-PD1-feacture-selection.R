library(magrittr) # pipe
library(DESeq2) # rse Ranged SummarizedExperiments
library(mlr) # machine learning wrapper
library(caret) # machine learning solution packages
library(glmnet) # elastic net and lasso general linear model
library(doParallel) # multiple parallel run elastic net
library(mRMRe) # mutual information
library(PCAtools) # PCA plot datasets
library(VennDiagram) # veendiagram
library(propagate) # big matrix correlation
library(biglasso) # big matrix lasso and elastic net

#filter gastric cancer,RNA-seq,anti-PD1,pretreatment metadata -------------
readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response) ->metadata
data.frame(Response=metadata$Response,batch=metadata$SRA_Study,row.names = metadata$Run)->Condition
Condition$Response %>%
  gsub("CR","R",.)%>%
  gsub("PR","R",.)%>%
  gsub("SD","NR",.)%>%
  gsub("PD","NR",.)->Condition$Response
Condition$batch %>%
  gsub("SRP070710",1,.)%>%
  gsub("SRP150548",2,.)%>%
  gsub("SRP094781",3,.)->Condition$batch

#make difference----------------------------
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1_pretreatment_Symbol_count_expr.txt",
           sep="\t",header = T,as.is = TRUE) ->melanoma_PD1_count_expr
melanoma.se <- SummarizedExperiment(assays = as.matrix(melanoma_PD1_count_expr), 
                                 colData = Condition)#,design= ~ Response + batch
# readr::write_rds(x = melanoma.se, path = '/data/liull/immune-checkpoint-blockade/melanoma_PD1.se.rds.gz', compress = 'gz')

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
  write.table(x = .mat, file = file.path('/data/liull/immune-checkpoint-blockade/machine_learning/', .filename), sep = '\t')
  if (file.exists(.filename)) file.remove(.filename)
  file.symlink(from = file.path('data', .filename), to = .filename)
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
    seed = 1234, ncores = 10, nfolds = 10
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

melanoma.feats <- fn_select_features(.se.name = "melanoma.se")

# doparallel --------------------------------------------------------------

cl <- makePSOCKcluster(names = 50)
registerDoParallel(cl = cl)
stopCluster(cl)
registerDoSEQ()
getDoParWorkers()


# function ----------------------------------------------------------------

# EDA
fn_pca_biplot_individual <- function(.se, .title = 'OC521') {
  .se.m <- assay(.se)
  .metadata <- colData(.se)[, c('oc', 'class'), drop = FALSE]
  .p <- pca(mat = .se.m, metadata = .metadata, removeVar = 0.2)
  .biplot <- biplot(
    .p,
    colby = 'class',
    lab = FALSE,
    legendPosition = 'right',
    title = glue::glue("{.title} PCA")
  )
  .biplot
}
fn_draw_venn <- function(.list){
  if (!is.null(dev.list())) dev.off()
  enn.plot <- venn.diagram(
    x = .list,
    filename = NULL,
    main = 'Feature Intersection',
    col = "transparent",
    fill = ggthemes::gdocs_pal()(length(.list)),
    alpha = 0.5,
    cex = 2.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.default.pos = "text",
    cat.cex = 2.5,
    cat.fontfamily = "serif",
    cat.pos = 0
  )
}



# For feature selection

 
fn_sd_median <- function(.df) {
  .sd <- apply(
    X = .df[, -1],
    MARGIN = 2,
    FUN = sd
  )
  .median <- apply(
    X = .df[, -1],
    MARGIN = 2,
    FUN = median
  )
  intersect(
    names(.sd[.sd > 0]),
    names(.median[.median > 0])
  ) -> .feats.sd.median
  .feats.sd.median
}




# Train and predict
fn_resample <- function(.se, .id) {
  .d <- cbind(t(assay(.se)), as.data.frame(colData(.se)[, 'class', drop = FALSE]))
  .task <- makeClassifTask(
    id = .id, data = .d,
    target = 'class', positive = 'M'
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
fn_caret_train_model <- function(.se, method = 'svmRadial') {
  .d <- fn_se2data.frame(.se = .se)
  .model <- caret::train(
    class ~ .,
    data = .d,
    method = method,
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
  .t <- as.numeric(.c)
  aList <- fn_getROC_AUC(probs = .pred$M, true_Y = .t)
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
fn_ensembl_predict <- function(.model, .se, .title = 'OC521') {
  .d.pred <- fn_predict(.model = .model, .se = .se)
  .d.pred.alist <- fn_alist(.c = .se$class, .pred = .d.pred)
  .d.p <- fn_plot_auc(.alist = .d.pred.alist, .t = .title)
  list(predict = .d.pred, alist = .d.pred.alist, plot = .d.p)
}

# feature selection -------------------------------------------------------

total816.feats <- fn_select_features(.se.name = "total816.se")
oc521.feats <- fn_select_features(.se.name = "oc521.se")
oc172.feats <- fn_select_features(.se.name = "oc172.se")
oc79.feats <- fn_select_features(.se.name = "oc79.se")
oc44.feats <- fn_select_features(.se.name = "oc44.se")
oc266.feats <- fn_select_features(.se.name = "oc266.se")
mrna.feats <- fn_select_features(.se.name = "mrna.se")
smart.feats <- fn_select_features(.se.name = "smart.se")

total816.feats <- readr::read_rds(path = 'data/total816.feats.rds.gz')
oc521.feats <- readr::read_rds(path = 'data/oc521.feats.rds.gz')
oc172.feats <- readr::read_rds(path = 'data/oc172.feats.rds.gz')
oc79.feats <- readr::read_rds(path = 'data/oc79.feats.rds.gz')
oc44.feats <- readr::read_rds(path = 'data/oc44.feats.rds.gz')
oc266.feats <- readr::read_rds(path = 'data/oc266.feats.rds.gz')
mrna.feats <- readr::read_rds(path = 'data/mrna.feats.rds.gz')
smart.feats <- readr::read_rds(path = 'data/smart.feats.rds.gz')

fn_resample(.se = total816.se[total816.feats, ], .id = 'total816')
fn_draw_venn(
  .list = list(
    total816 = total816.feats,
    oc521 = oc521.feats,
    lxc_ens = lxc_ens
  )
) %>% grid::grid.draw()

smart.model <- fn_caret_train_model(.se = oc521.se[total816.old.feats, ])
panel.smart.test.se.auc <- fn_ensembl_predict(.model = smart.model, .se = oc266.se[total816.old.feats, ], .title = 'smart.test')

panel.smart.test.se.auc$plot

fn_pca_biplot_individual(.se = mrna.se, .title = 'mRNA')
fn_pca_biplot_individual(.se = smart.se, .title = 'SMART')






