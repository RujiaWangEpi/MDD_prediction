# Models

Scripts for predicting major depressive disorder (MDD) diagnosis, number of MDD episodes and age of MDD onset in GLAD+ and UK Biobank (UKB). Predictors are family history, polygenic risk scores (PRS), childhood trauma and background factors (sex and ten genetic principal components).

## Scripts

| Script | What it does |
|---|---|
| `Model_builder_V2_nested_wt.R` | Nested cross-validated elastic net with class weights (5 random outer folds for evaluation; 10 inner folds, stratified by case status, for tuning). Saves the outer-fold models. |
| `Model_builder_V2_nested_wt_refit.R` | Same nested elastic net, plus a refit of each model on the full development sample. It also writes calibration, AUC/R² confidence intervals and PPV/NPV from the pooled out-of-fold predictions. |
| `External_validation_fold.R` | External validation with one outer-fold model (`--fold k`) from one cohort, applied unchanged to the other cohort. |
| `External_validation_refit.R` | External validation with the full-sample refit model from one cohort, applied unchanged to the other cohort. |
| `prediction_metrics.R` | Shared functions for AUC, liability-scale R², calibration, decile tables, PPV/NPV and calibration plots. Keep it in the same folder as the scripts above. |
| `Compare_AUC_DeLong.R` | Paired DeLong tests comparing AUCs between predictor groups, with FDR and Bonferroni correction. |
| `xgb_mdd.R` | XGBoost sensitivity analysis: 5-fold cross-validation within each cohort and external validation between cohorts, with calibration and SHAP values. |
| `Cox_mdd_inc.R` | Cox proportional hazards model for incident MDD in UKB, with a proportional-hazards check. |
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

# 3. Compare AUCs between predictor groups
Rscript Compare_AUC_DeLong.R

# 4. XGBoost sensitivity analysis
Rscript xgb_mdd.R --outdir results/xgb
```

## Notes

- `--outcome_pop_prev` is the population prevalence used for the liability-scale R² (5% for GLAD+, 15% for UKB in the paper).
- For external validation, the target cohort's predictors file must list the same predictor files, in the same order and with the same column names, as in the development cohort.
- Add `--keep ids.txt` to the external validation scripts to restrict the target cohort, for example to the resampled UKB sample.
- To apply the MDD diagnosis model to a severity outcome, give the number-of-episodes or age-of-onset file as `--pheno`.

## Requirements

R (4.3 or later) with `optparse`, `data.table`, `glmnet`, `caret`, `doMC`, `pROC`, `psych`, `verification` and `ggplot2`. `xgb_mdd.R` also needs `xgboost` and `Matrix`; `Cox_mdd_inc.R` needs `survival`.
