# Compare AUCs between predictor groups using paired DeLong tests (Table S9; former Table S7)
# Input: pooled out-of-fold ROC objects written by Model_builder_V2_nested_wt_refit.R:
#   <out>.<group>.Indep_auc.rds   (one per predictor group, plus "All")
# All groups are evaluated on identical outer folds (fixed seed, same analytic sample), so the
# ROC objects are paired. Comparisons with the full model are nested; following Demler et al.
# (Stat Med 2012) these are interpreted alongside the AUC difference (dAUC column).
# Usage: Rscript Compare_AUC_DeLong.R <glad_out_prefix> <ukb_out_prefix> <output_file>
library(pROC)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript Compare_AUC_DeLong.R <glad_out_prefix> <ukb_out_prefix> <output_file>")

# Group labels as used in the predictors file (group column) of the development runs
grp <- c(FH = "fh", PRS = "prs", CT = "ct", Demo = "demo", All = "All")

compare_auc <- function(prefix) {
  files <- setNames(paste0(prefix, ".", grp, ".Indep_auc.rds"), names(grp))
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) stop("ROC objects not found:\n  ", paste(missing, collapse = "\n  "))
  roc_list <- lapply(files, readRDS)
  pairs <- combn(names(roc_list), 2, simplify = FALSE)
  res <- do.call(rbind, lapply(pairs, function(x) {
    stopifnot(identical(roc_list[[x[1]]]$response, roc_list[[x[2]]]$response))  # paired on the same participants
    a1 <- as.numeric(auc(roc_list[[x[1]]])); a2 <- as.numeric(auc(roc_list[[x[2]]]))
    data.frame(model1 = x[1], model2 = x[2], auc1 = a1, auc2 = a2, dAUC = a2 - a1,
               p_value = roc.test(roc_list[[x[1]]], roc_list[[x[2]]],
                                  method = "delong", paired = TRUE)$p.value)
  }))
  res$p_fdr  <- p.adjust(res$p_value, method = "fdr")
  res$p_bonf <- p.adjust(res$p_value, method = "bonferroni")
  res
}

glad <- compare_auc(args[1])
ukb  <- compare_auc(args[2])
write.table(rbind(cbind(Cohort = "GLAD+", glad), cbind(Cohort = "UKB", ukb)),
            args[3], sep = "\t", row.names = FALSE, quote = FALSE)
