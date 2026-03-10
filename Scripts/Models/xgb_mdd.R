#### XGBoost for MDD prediction in GLAD+ and UKB ####
#### Rujia Wang 2026-01-30 ####

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(xgboost)
  library(pROC)
})

# -------------------------
# Args (no extra packages)
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

out_dir <- get_arg("--outdir", default = "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ID_COL <- "IID"
Y_COL  <- "Outcome_var"
K_FOLDS <- 5
SEED <- 2026

`%||%` <- function(a, b) if (!is.null(a)) a else b

get_nthread <- function(default = 8L) {
  cand <- c("SLURM_CPUS_PER_TASK", "DX_CPUS", "DX_JOB_CPUS", "OMP_NUM_THREADS")
  for (nm in cand) {
    v <- suppressWarnings(as.integer(Sys.getenv(nm, NA)))
    if (!is.na(v) && v >= 1L) return(v)
  }
  default
}

prep_binary_y <- function(y) {
  if (is.factor(y)) y <- as.character(y)
  if (is.character(y)) {
    y_low <- tolower(y)
    if (all(y_low %in% c("0","1","case","control","yes","no","true","false"))) {
      y <- fifelse(y_low %in% c("1","case","yes","true"), 1L, 0L)
    } else stop("Outcome_var is character with unexpected levels. Recode to 0/1.")
  } else {
    y <- as.integer(y)
    if (all(na.omit(unique(y)) %in% c(1L, 2L))) y <- y - 1L
  }
  stopifnot(all(na.omit(unique(y)) %in% c(0L, 1L)))
  as.integer(y)
}

make_stratified_folds <- function(y, k = 5, seed = 1) {
  set.seed(seed)
  y <- as.integer(y)
  idx1 <- sample(which(y == 1L))
  idx0 <- sample(which(y == 0L))
  folds <- vector("list", k)
  for (i in seq_len(k)) {
    folds[[i]] <- c(idx1[seq(i, length(idx1), by = k)],
                    idx0[seq(i, length(idx0), by = k)])
  }
  folds
}

make_sparse_numeric <- function(df, y_col, id_col, feature_cols = NULL) {
  df <- as.data.table(df)
  y <- as.numeric(df[[y_col]])
  
  if (is.null(feature_cols)) feature_cols <- setdiff(names(df), c(y_col, id_col))
  Xdf <- df[, ..feature_cols]
  for (cc in feature_cols) Xdf[[cc]] <- as.numeric(Xdf[[cc]])
  
  X <- Matrix(as.matrix(Xdf), sparse = TRUE)
  colnames(X) <- feature_cols
  list(X = X, y = y, feature_cols = feature_cols)
}

calc_auc <- function(y_true, y_prob) {
  as.numeric(pROC::auc(pROC::roc(y_true, y_prob, quiet = TRUE)))
}

calc_calibration <- function(y_true, p_hat) {
  eps <- 1e-6
  p_hat <- pmin(pmax(p_hat, eps), 1 - eps)
  fit <- glm(y_true ~ qlogis(p_hat), family = binomial())
  data.table(
    cal_intercept = unname(coef(fit)[1]),
    cal_slope     = unname(coef(fit)[2]),
    brier         = mean((p_hat - y_true)^2)
  )
}

# -------------------------
# Save ROC + AUC as RDS (+ small CSV)
# -------------------------
save_roc_auc_rds <- function(y_true, p_hat, prefix, meta = list()) {
  y_true <- prep_binary_y(y_true)
  roc_obj <- pROC::roc(response = y_true, predictor = p_hat, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  out <- c(list(roc = roc_obj, auc = auc_val, y_true = y_true, p_hat = p_hat), meta)
  saveRDS(out, file.path(out_dir, paste0(prefix, "_roc_auc.rds")))
  
  fwrite(data.table(prefix = prefix, auc = auc_val),
         file.path(out_dir, paste0(prefix, "_auc.csv")))
  invisible(out)
}

# -------------------------
# SHAP summary (ALL features; summary only)
# -------------------------
xgb_shap_summary_all <- function(model, dmat, feature_names = NULL,
                                 include_quantiles = TRUE,
                                 probs = c(0.01, 0.05, 0.50, 0.95, 0.99)) {
  shap_mat <- predict(model, dmat, predcontrib = TRUE)
  shap_df <- as.data.frame(shap_mat)
  if ("BIAS" %in% names(shap_df)) shap_df$BIAS <- NULL
  if (!is.null(feature_names) && length(feature_names) == ncol(shap_df)) {
    names(shap_df) <- feature_names
  }
  feats <- names(shap_df)
  
  out <- rbindlist(lapply(feats, function(f) {
    v <- shap_df[[f]]
    dt <- data.table(
      feature = f,
      mean_abs_shap = mean(abs(v), na.rm = TRUE),
      mean_shap     = mean(v, na.rm = TRUE),
      sd_shap       = sd(v, na.rm = TRUE)
    )
    if (include_quantiles) {
      qs <- as.numeric(stats::quantile(v, probs = probs, na.rm = TRUE))
      qnames <- paste0("q", sprintf("%02d", as.integer(probs * 100)))
      dt[, (qnames) := as.list(qs)]
    }
    dt
  }))
  setorder(out, -mean_abs_shap)
  out
}

# -------------------------
# Tuning grid
# -------------------------
param_grid <- expand.grid(
  eta              = c(0.03, 0.05, 0.1),
  max_depth        = c(3, 4, 6),
  min_child_weight = c(1, 5),
  subsample        = 0.8,
  colsample_bytree = 0.8
)

# -------------------------
# Tune with xgb.cv (robust fallback)
# -------------------------
tune_xgb_cv <- function(X, y, folds, grid,
                        nrounds_max = 3000,
                        early_stopping_rounds = 50,
                        seed = 1) {
  
  pos <- sum(y == 1); neg <- sum(y == 0)
  spw <- if (pos > 0) neg / pos else 1
  nthread <- get_nthread(8L)
  d <- xgb.DMatrix(X, label = y, missing = NA)
  
  best_score <- -Inf
  best_params <- NULL
  best_nrounds <- NULL
  
  for (i in seq_len(nrow(grid))) {
    p <- list(
      booster = "gbtree",
      objective = "binary:logistic",
      eval_metric = "auc",
      nthread = nthread,
      eta = grid$eta[i],
      max_depth = grid$max_depth[i],
      min_child_weight = grid$min_child_weight[i],
      subsample = grid$subsample[i],
      colsample_bytree = grid$colsample_bytree[i],
      gamma = 0,
      lambda = 1,
      alpha = 0,
      scale_pos_weight = spw
    )
    
    set.seed(seed + i)
    cv <- tryCatch(
      xgb.cv(
        params = p,
        data = d,
        folds = folds,
        nrounds = nrounds_max,
        early_stopping_rounds = early_stopping_rounds,
        maximize = TRUE,
        verbose = 0
      ),
      error = function(e) NULL
    )
    
    if (is.null(cv) || is.null(cv$evaluation_log) || nrow(cv$evaluation_log) == 0) next
    if (is.null(cv$best_iteration) || length(cv$best_iteration) == 0) next
    
    score <- max(cv$evaluation_log$test_auc_mean)
    nbest <- cv$best_iteration
    
    if (score > best_score) {
      best_score <- score
      best_params <- p
      best_nrounds <- nbest
    }
  }
  
  if (is.null(best_nrounds) || length(best_nrounds) == 0) {
    message("WARNING: tuning failed; using default params and nrounds=200.")
    best_params <- list(
      booster="gbtree", objective="binary:logistic", eval_metric="auc",
      nthread=nthread, eta=0.05, max_depth=4, min_child_weight=1,
      subsample=0.8, colsample_bytree=0.8, gamma=0, lambda=1, alpha=0,
      scale_pos_weight=spw
    )
    best_nrounds <- 200L
    best_score <- NA_real_
  }
  
  list(best_params = best_params, best_nrounds = best_nrounds, tuning_cv_mean_auc = best_score)
}

# ============================================================
# 5-fold CV (OOF pred) + calibration + SHAP + ROC/AUC saving
# ============================================================
xgb_cv_5fold <- function(df, dataset_label, k = 5, seed = 2026, do_shap = TRUE) {
  df <- as.data.table(df)
  df[[Y_COL]] <- prep_binary_y(df[[Y_COL]])
  df <- df[is.finite(df[[Y_COL]])]
  if (length(unique(df[[Y_COL]])) < 2) stop(dataset_label, ": outcome has <2 classes.")
  
  feature_cols <- setdiff(names(df), c(Y_COL, ID_COL))
  mm <- make_sparse_numeric(df, Y_COL, ID_COL, feature_cols)
  X <- mm$X
  y <- as.integer(mm$y)
  
  folds <- make_stratified_folds(y, k = k, seed = seed)
  tuned <- tune_xgb_cv(X, y, folds, param_grid, seed = seed)
  
  fold_auc <- numeric(k)
  oof_pred <- rep(NA_real_, length(y))
  
  for (j in seq_len(k)) {
    test_idx  <- folds[[j]]
    train_idx <- setdiff(seq_len(nrow(df)), test_idx)
    
    dtrain <- xgb.DMatrix(X[train_idx, , drop = FALSE], label = y[train_idx], missing = NA)
    dtest  <- xgb.DMatrix(X[test_idx,  , drop = FALSE], missing = NA)
    
    model <- xgb.train(params = tuned$best_params, data = dtrain, nrounds = tuned$best_nrounds, verbose = 0)
    p_hat <- predict(model, dtest)
    
    oof_pred[test_idx] <- p_hat
    fold_auc[j] <- calc_auc(y[test_idx], p_hat)
    message(sprintf("[%s] Fold %d/%d AUC=%.4f", dataset_label, j, k, fold_auc[j]))
  }
  
  cal_oof <- calc_calibration(y, oof_pred)
  
  final_model <- xgb.train(
    params = tuned$best_params,
    data = xgb.DMatrix(X, label = y, missing = NA),
    nrounds = tuned$best_nrounds,
    verbose = 0
  )
  
  shap_all <- NULL
  if (isTRUE(do_shap)) {
    shap_all <- xgb_shap_summary_all(
      final_model,
      xgb.DMatrix(X, label = y, missing = NA),
      feature_names = feature_cols,
      include_quantiles = TRUE
    )
    shap_all[, `:=`(dataset = dataset_label, model = "XGBoost", source = "CV_full_train_model")]
  }
  
  list(
    dataset = dataset_label,
    model = "XGBoost",
    cv_folds = k,
    seed = seed,
    nthread = tuned$best_params$nthread %||% get_nthread(8L),
    best_params = tuned$best_params,
    best_nrounds = tuned$best_nrounds,
    tuning_cv_mean_auc = tuned$tuning_cv_mean_auc,
    fold_auc = fold_auc,
    auc_mean = mean(fold_auc),
    auc_sd = sd(fold_auc),
    y_true = y,
    oof_pred = oof_pred,
    cal_oof = cbind(data.table(dataset = dataset_label, model = "XGBoost", source = "OOF_CV"), cal_oof),
    feature_names = feature_cols,
    shap_all = shap_all
  )
}

# ============================================================
# External validation + calibration + SHAP + ROC/AUC saving
# ============================================================
xgb_external <- function(train_df, test_df, train_label, test_label, k = 5, seed = 2026, do_shap = TRUE) {
  train_df <- as.data.table(train_df)
  test_df  <- as.data.table(test_df)
  
  train_df[[Y_COL]] <- prep_binary_y(train_df[[Y_COL]])
  test_df[[Y_COL]]  <- prep_binary_y(test_df[[Y_COL]])
  
  train_df <- train_df[is.finite(train_df[[Y_COL]])]
  test_df  <- test_df[is.finite(test_df[[Y_COL]])]
  
  feature_cols <- setdiff(names(train_df), c(Y_COL, ID_COL))
  miss <- setdiff(feature_cols, names(test_df))
  if (length(miss) > 0) for (cc in miss) test_df[[cc]] <- NA_real_
  
  mm_tr <- make_sparse_numeric(train_df, Y_COL, ID_COL, feature_cols)
  mm_te <- make_sparse_numeric(test_df,  Y_COL, ID_COL, feature_cols)
  
  Xtr <- mm_tr$X; ytr <- as.integer(mm_tr$y)
  Xte <- mm_te$X; yte <- as.integer(mm_te$y)
  
  folds <- make_stratified_folds(ytr, k = k, seed = seed)
  tuned <- tune_xgb_cv(Xtr, ytr, folds, param_grid, seed = seed)
  
  model <- xgb.train(
    params = tuned$best_params,
    data = xgb.DMatrix(Xtr, label = ytr, missing = NA),
    nrounds = tuned$best_nrounds,
    verbose = 0
  )
  
  p_hat <- predict(model, xgb.DMatrix(Xte, missing = NA))
  auc_ext <- calc_auc(yte, p_hat)
  cal_ext <- calc_calibration(yte, p_hat)
  
  shap_all <- NULL
  if (isTRUE(do_shap)) {
    shap_all <- xgb_shap_summary_all(
      model,
      xgb.DMatrix(Xtr, label = ytr, missing = NA),
      feature_names = feature_cols,
      include_quantiles = TRUE
    )
    shap_all[, `:=`(train_dataset = train_label, test_dataset = test_label,
                    model = "XGBoost", source = "External_model_SHAP_on_train")]
  }
  
  list(
    train_dataset = train_label,
    test_dataset  = test_label,
    model = "XGBoost",
    cv_folds_train = k,
    seed = seed,
    nthread = tuned$best_params$nthread %||% get_nthread(8L),
    best_params = tuned$best_params,
    best_nrounds = tuned$best_nrounds,
    tuning_cv_mean_auc_train = tuned$tuning_cv_mean_auc,
    external_auc = auc_ext,
    y_true = yte,
    pred = p_hat,
    cal_external = cbind(data.table(train_dataset = train_label, test_dataset = test_label,
                                    model = "XGBoost", source = "External"),
                         cal_ext),
    shap_all = shap_all
  )
}

# -------------------------
# Save outputs (CV)
# -------------------------
save_cv_outputs <- function(obj, prefix) {
  saveRDS(obj, file.path(out_dir, paste0(prefix, ".rds")))
  
  fwrite(data.table(dataset=obj$dataset, model=obj$model, cv_folds=obj$cv_folds,
                    auc_mean=obj$auc_mean, auc_sd=obj$auc_sd,
                    tuning_cv_mean_auc=obj$tuning_cv_mean_auc,
                    best_nrounds=obj$best_nrounds, nthread=obj$nthread),
         file.path(out_dir, paste0(prefix, "_cv_summary.csv")))
  
  fwrite(data.table(dataset=obj$dataset, model=obj$model, fold=seq_along(obj$fold_auc), auc=obj$fold_auc),
         file.path(out_dir, paste0(prefix, "_cv_fold_auc.csv")))
  
  fwrite(obj$cal_oof, file.path(out_dir, paste0(prefix, "_cv_calibration_oof.csv")))
  
  params_dt <- as.data.table(obj$best_params)
  params_dt[, `:=`(dataset=obj$dataset, model=obj$model, best_nrounds=obj$best_nrounds)]
  fwrite(params_dt, file.path(out_dir, paste0(prefix, "_cv_best_params.csv")))
  
  if (!is.null(obj$shap_all) && nrow(obj$shap_all) > 0) {
    fwrite(obj$shap_all, file.path(out_dir, paste0(prefix, "_shap_all.csv")))
  }
  
  # ROC/AUC RDS using OOF predictions
  save_roc_auc_rds(
    y_true = obj$y_true,
    p_hat  = obj$oof_pred,
    prefix = paste0(prefix, "_oof"),
    meta = list(dataset = obj$dataset, model = obj$model, source = "OOF_CV")
  )
}

# -------------------------
# Save outputs (External)
# -------------------------
save_external_outputs <- function(obj, prefix) {
  saveRDS(obj, file.path(out_dir, paste0(prefix, ".rds")))
  
  fwrite(data.table(train_dataset=obj$train_dataset, test_dataset=obj$test_dataset,
                    model=obj$model, cv_folds_train=obj$cv_folds_train,
                    external_auc=obj$external_auc,
                    tuning_cv_mean_auc_train=obj$tuning_cv_mean_auc_train,
                    best_nrounds=obj$best_nrounds, nthread=obj$nthread),
         file.path(out_dir, paste0(prefix, "_external_summary.csv")))
  
  fwrite(obj$cal_external, file.path(out_dir, paste0(prefix, "_external_calibration.csv")))
  
  params_dt <- as.data.table(obj$best_params)
  params_dt[, `:=`(train_dataset=obj$train_dataset, test_dataset=obj$test_dataset,
                   model=obj$model, best_nrounds=obj$best_nrounds)]
  fwrite(params_dt, file.path(out_dir, paste0(prefix, "_external_best_params.csv")))
  
  if (!is.null(obj$shap_all) && nrow(obj$shap_all) > 0) {
    fwrite(obj$shap_all, file.path(out_dir, paste0(prefix, "_shap_all.csv")))
  }
  
  # ROC/AUC RDS using external predictions
  save_roc_auc_rds(
    y_true = obj$y_true,
    p_hat  = obj$pred,
    prefix = paste0(prefix, "_external"),
    meta = list(train_dataset = obj$train_dataset, test_dataset = obj$test_dataset,
                model = obj$model, source = "External")
  )
}

# -------------------------
# Load inputs (your current approach)
# -------------------------
load("GLAD_replicated.Rdata")          
load("UKB_replicated.Rdata")           
load("UKB_replicated_resample.Rdata")  

mdd_onset<-readRDS("mdd_onset_cox.rds")
mdd_onset$Outcome_var<-mdd_onset$obs
mdd_onset1<-mdd_onset[,c(2,56,6:33)]

# -------------------------
# Run: within-dataset CV
# -------------------------
glad_cv <- xgb_cv_5fold(GLAD_test, "GLAD+", k = K_FOLDS, seed = SEED, do_shap = TRUE)
save_cv_outputs(glad_cv, "glad_xgb")

ukb_cv <- xgb_cv_5fold(UKB_test, "UKB", k = K_FOLDS, seed = SEED, do_shap = TRUE)
save_cv_outputs(ukb_cv, "ukb_xgb")

ukb_inc_cv <- xgb_cv_5fold(mdd_onset1, "UKB_inc", k = K_FOLDS, seed = SEED, do_shap = TRUE)
save_cv_outputs(ukb_inc_cv, "ukb_inc_xgb")

# -------------------------
# Run: external validation
# -------------------------
glad_to_ukb <- xgb_external(GLAD_test, UKB_test, "GLAD+", "UKB", k = K_FOLDS, seed = SEED, do_shap = TRUE)
save_external_outputs(glad_to_ukb, "glad_to_ukb_xgb")

glad_to_ukb_resample <- xgb_external(GLAD_test, UKB_test_resample, "GLAD+", "UKB_resample", k = K_FOLDS, seed = SEED, do_shap = TRUE)
save_external_outputs(glad_to_ukb_resample, "glad_to_ukb_resample_xgb")

ukb_to_glad <- xgb_external(UKB_test, GLAD_test, "UKB", "GLAD+", k = K_FOLDS, seed = SEED, do_shap = TRUE)
save_external_outputs(ukb_to_glad, "ukb_to_glad_xgb")

message("DONE. Outputs in: ", normalizePath(out_dir))

q(save = "no")
