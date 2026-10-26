# Compare AUCs between predictor groups using paired DeLong tests (Table S7)
# Input: pooled out-of-fold ROC objects (*.Indep_auc_<group>.rds) from Model_builder_V2_nested_wt.R
library(pROC)

compare_auc <- function(files) {
  roc_list <- lapply(files, readRDS)
  pairs <- combn(names(roc_list), 2, simplify = FALSE)
  res <- do.call(rbind, lapply(pairs, function(x) {
    data.frame(model1 = x[1], model2 = x[2],
               auc1 = as.numeric(auc(roc_list[[x[1]]])),
               auc2 = as.numeric(auc(roc_list[[x[2]]])),
               p_value = roc.test(roc_list[[x[1]]], roc_list[[x[2]]],
                                  method = "delong", paired = TRUE)$p.value)
  }))
  res$p_fdr  <- p.adjust(res$p_value, method = "fdr")
  res$p_bonf <- p.adjust(res$p_value, method = "bonferroni")
  res
}

grp <- c(FH = "fh", PRS = "prs", CT = "ct", Demo = "demo", All = "all")
glad <- compare_auc(setNames(paste0("Results/GLAD/mdd_all_pred.Indep_auc_", grp, ".rds"), names(grp)))
ukb  <- compare_auc(setNames(paste0("Results/UKB/nested_ukb_mdd_allpred.Indep_auc_", grp, ".rds"), names(grp)))

write.table(rbind(cbind(Cohort = "GLAD+", glad), cbind(Cohort = "UKB", ukb)),
            "Results/TableS7_DeLong.txt", sep = "\t", row.names = FALSE, quote = FALSE)