#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --output=Nested_cv_mdd_ukb.log
#SBATCH --export=NONE
#SBATCH --get-user-env=60L
#SBATCH --begin=now+0hour

module load r/4.3.0-gcc-13.2.0-withx-rmath-standalone-python-3.11.6

cd /path/to/project   # edit

# MDD diagnosis

# --assoc F: the per-predictor univariate associations are a legacy pipeline output and are not used for any reported result.
# --compare_predictors F: single-predictor models are not reported (predictor contributions are taken from the full model; Figure S2).
Rscript Scripts/Models/Model_builder_V2_nested_wt_refit.R \
    --pheno predictors/UKB_MDD.txt  \
    --out results/ukb_mdd/nested_ukb_mdd_allpred \
    --n_core 1 \
    --compare_predictors F \
    --assoc F \
    --outcome_pop_prev 0.15 \
    --predictors predictors/predictors_map_ukb.txt

# MDD episodes 

Rscript Scripts/Models/Model_builder_V2_nested_wt_refit.R \
    --pheno predictors/UKB_MDD_episodes.txt  \
    --out results/ukb_episodes/nested_ukb_mdd_episodes \
    --n_core 1 \
    --compare_predictors F \
    --assoc F \
    --predictors predictors/predictors_map_ukb.txt

