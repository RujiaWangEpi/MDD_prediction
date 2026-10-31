# Models

Scripts for predicting major depressive disorder (MDD) diagnosis, number of MDD episodes and age of MDD onset in GLAD+ and UK Biobank (UKB). Predictors are family history, polygenic risk scores (PRS), childhood trauma and background factors (sex and ten genetic principal components).

## Scripts

| Script | What it does |
|---|---|
| `Model_builder_V2_nested_wt.R` | Nested cross-validated elastic net with class weights (5 random outer folds for evaluation; 10 inner folds, stratified by case status, for tuning). Saves the outer-fold models. |
| `Model_builder_V2_nested_wt_refit.R` | Same nested elastic net, plus a refit of each model on the full development sample. Writes outer-fold and refit coefficients (Tables S7, S7a; Figures S2-S4), fold AUCs and selected tuning parameters (Tables S6a, S9-S13). Outputs are named by predictor group. |
| `External_validation_fold.R` | External validation with one outer-fold model (`--fold k`) from one cohort, applied unchanged to the other cohort. |
| `External_validation_refit.R` | External validation with the full-sample refit model from one cohort, applied unchanged to the other cohort. |
| `prediction_metrics.R` | Shared functions for AUC (Hanley-McNeil and DeLong CIs), observed and liability-scale R² (Table S8 across K), calibration (Brier, slope, intercept), decile tables, PPV/NPV at cut-offs 0.3-0.7 and expected PPV/NPV by prevalence (Tables S14-S18, S20). Used by the external validation scripts; apply it to `<out>.<group>.indep_pred.rds` for the within-cohort rows. |
| `Compare_AUC_DeLong.R` | Paired DeLong tests comparing AUCs between predictor groups on the pooled out-of-fold ROC objects, with ΔAUC, FDR and Bonferroni correction (Table S9). |
| `xgb_mdd.R` | XGBoost sensitivity analysis: 5-fold stratified cross-validation within each cohort and external validation between cohorts, with calibration (Figure 3C-F, Table S14). In the reported analyses the grid search did not return a result and prespecified hyperparameters were used (see the note at the top of the script and Table S6); `--no_tune` reproduces this. |
| `Nested_prediction.sh` | Example SLURM job for running the nested elastic net on an HPC cluster. |

## Typical workflow

```bash
# 1. Build and internally validate the models in each cohort
Rscript Model_builder_V2_nested_wt_refit.R \
    --pheno predictors/GLAD_MDD.txt \
    --predictors predictors/predictors_map_glad.txt \
    --outcome_pop_prev 0.05 \
    --out results/glad_mdd/glad_mdd

# 2. External validation, with either an outer-fold model ...
Rscript External_validation_fold.R \
    --model_prefix results/glad_mdd/glad_mdd --fold 1 \
    --pheno predictors/UKB_MDD.txt \
    --predictors predictors/predictors_map_ukb.txt \
    --outcome_pop_prev 0.15 \
    --train_label GLAD+ --test_label UKB \
    --out results/external/glad_to_ukb_fold1

# ... or the full-sample refit model
Rscript External_validation_refit.R \
    --model_prefix results/glad_mdd/glad_mdd \
    --pheno predictors/UKB_MDD.txt \
    --predictors predictors/predictors_map_ukb.txt \
    --outcome_pop_prev 0.15 \
    --train_label GLAD+ --test_label UKB \
    --out results/external/glad_to_ukb_refit

# 3. Compare AUCs between predictor groups (Table S9)
Rscript Compare_AUC_DeLong.R results/glad_mdd/glad_mdd results/ukb_mdd/ukb_mdd results/TableS9_DeLong.txt

# 4. XGBoost sensitivity analysis (as reported: prespecified hyperparameters)
Rscript xgb_mdd.R --outdir results/xgb --no_tune
```

## Analyses and where they are produced

| Manuscript item | Script / step |
|---|---|
| Table S9 (group AUCs, DeLong), Table S10 (MDD PRS vs mPRS), Table S11 (without mPRS), Table S12 (incremental over background), Table S13 (22 traits) | `Model_builder_V2_nested_wt_refit.R` with the corresponding predictors file (one row per predictor file, `group` column); `Compare_AUC_DeLong.R` |
| Tables S7, S7a; Figures S2-S4 (coefficients, fold selection) | `<out>.<group>.coefficients.txt` from `Model_builder_V2_nested_wt_refit.R` (standardised coefficient = coefficient x predictor SD) |
| Table S8 (liability R² across K), Tables S14-S18, S20 | `External_validation_fold.R` / `External_validation_refit.R` (external rows); `prediction_metrics.R` on pooled out-of-fold predictions (within-cohort rows) |
| Resampled UKB (Table S20): cases sampled to match the GLAD+ episode distribution plus 845 random controls, seed 2026 | `--keep` file for the external validation scripts (ID list; the sampling itself is a few lines of `sample()` by episode category and is documented in the Supplementary Methods) |
| Figure 3C-F, Table S14 XGBoost rows | `xgb_mdd.R` |
| Figures 2 and 4, Figure S5-S6 | plotted from `<out>.ext_deciles.txt` / `<out>.ext_predictions.txt` with ggplot2 (plotting scripts not included) |

## Notes

- `--outcome_pop_prev` is the population prevalence used for the liability-scale R² (5% for GLAD+, 15% for UKB in the paper).
- For external validation, the target cohort's predictors file must list the same predictor files, in the same order and with the same column names, as in the development cohort.
- Add `--keep ids.txt` to the external validation scripts to restrict the target cohort, for example to the resampled UKB sample.
- To apply the MDD diagnosis model to a severity outcome, give the number-of-episodes or age-of-onset file as `--pheno`.

## Requirements

R (4.3 or later) with `optparse`, `data.table`, `glmnet`, `caret`, `doMC`, `pROC`, `psych`, `verification` and `ggplot2`. `xgb_mdd.R` also needs `xgboost` and `Matrix`. Package versions used for the reported analyses are listed in Supplementary Table S21.
