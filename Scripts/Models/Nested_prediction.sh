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

cd project folder

# MDD diagnosis

Rscript scripts/Model_builder_V2_nested_wt.R \
    --pheno predictors/MDD_rm_episodes.txt  \
    --out results/mdd_wt_20250120/nested_ukb_mdd_allpred \
    --n_core 1 \
    --compare_predictors T \
    --assoc T \
    --outcome_pop_prev 0.15 \
    --predictors /predictors/predictors_map_gald.txt

# MDD episodes 

Rscript scripts/Model_builder_V2_nested.R \
    --pheno predictors/MDD_episodes.txt  \
    --out results/mdd_episodes_20250120/ukb_mdd_episode \
    --n_core 1 \
    --compare_predictors T \
    --assoc T \
    --predictors predictors/predictors_map_gald.txt

