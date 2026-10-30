#!/usr/bin/Rscript
# External_validation_fold.R
# Rujia Wang, 2026
#
# External validation using ONE outer-fold model: the elastic net model trained on 80% of the
# development cohort in the chosen outer fold of the nested cross-validation
# (Model_builder_V2_nested_wt_refit.R, or the earlier Model_builder_V2_nested_wt.R),
# applied unchanged to the other cohort (GLAD+ -> UKB, UKB -> GLAD+, GLAD+ -> resampled UKB).
#
# Outputs (prefix --out):
#   .ext_eval.txt             AUC (Hanley-McNeil and DeLong 95% CI), observed and liability R2 (95% CI),
#                             Brier score, calibration slope and intercept                         (Tables S14, S20)
#                             For a continuous outcome: standardised beta, observed R2 (95% CI)     (Tables S15, S20)
#   .ext_thresholds.txt       PPV, NPV, sensitivity, specificity at cut-offs 0.3-0.7                (Table S16)
#   .ext_expected_ppv_npv.txt expected PPV/NPV at other prevalences, from Se/Sp at 0.5             (Table S17)
#   .ext_deciles.txt          observed vs predicted by decile of predicted risk                     (Figure 3E-F, Figure S5)
#   .ext_calibration.pdf      calibration plot
#   .ext_roc.rds              pROC object (for Figure 3C-D and DeLong tests)
#   .ext_predictions.txt      individual-level predictions
#
# The target cohort must be prepared exactly like the development cohort: same predictors file
# (same files in the same order, same column names and coding, PRSs standardised to the same reference).
# Evaluation functions are in prediction_metrics.R, which must sit in the same folder as this script.
#
# Example (GLAD+ outer-fold 1 model applied to UKB, MDD diagnosis):
#   Rscript External_validation_fold.R \
#     --model_prefix results/glad_mdd/nested_glad_mdd_allpred \
#     --fold 1 \
#     --pheno predictors/UKB_MDD.txt \
#     --predictors predictors/predictors_map_ukb.txt \
#     --outcome_pop_prev 0.15 \
#     --train_label GLAD+ --test_label UKB \
#     --out results/external/glad_to_ukb_mdd_fold1

start.time <- Sys.time()
suppressMessages(library(optparse))

option_list = list(
make_option("--model_prefix", action="store", default=NA, type='character',
		help="--out prefix used when building the models in the development cohort [required]"),
make_option("--group", action="store", default='All', type='character',
		help="Model (predictor group) to validate, e.g. All, fh, megaPRS [default All]"),
make_option("--pheno", action="store", default=NA, type='character',
		help="Phenotype file for the target cohort (FID, IID, outcome) [required]"),
make_option("--predictors", action="store", default=NA, type='character',
		help="Predictors file for the target cohort, same format and order as in the development cohort [required]"),
make_option("--compare_predictors", action="store", default=F, type='logical',
		help="Use the same value as in the development run (affects predictor names) [optional]"),
make_option("--keep", action="store", default=NA, type='character',
		help="File listing individuals (FID, IID) to include, e.g. the resampled UKB subset [optional]"),
make_option("--fold", action="store", default=1, type='numeric',
		help="Outer fold whose model is applied to the target cohort [default 1]"),
make_option("--outcome_pop_prev", action="store", default=NA, type='numeric',
		help="Population prevalence for the liability-scale R2 in the target cohort (e.g. 0.15 UKB, 0.05 GLAD+) [optional]"),
make_option("--ppv_prev", action="store", default='0.10,0.20,0.36,0.69', type='character',
		help="Prevalences for expected PPV/NPV at the 0.5 cut-off [default 0.10,0.20,0.36,0.69]"),
make_option("--train_label", action="store", default='Development', type='character',
		help="Label for the development cohort [optional]"),
make_option("--test_label", action="store", default='Target', type='character',
		help="Label for the target cohort [optional]"),
make_option("--out", action="store", default=NA, type='character',
		help="Prefix for output files [required]")
)
opt = parse_args(OptionParser(option_list=option_list))
for(req in c('model_prefix','pheno','predictors','out')) if(is.na(opt[[req]])) stop('--', req, ' is required')
dir.create(dirname(opt$out), showWarnings=FALSE, recursive=TRUE)

suppressMessages(library(data.table))
suppressMessages(library(glmnet))
suppressMessages(library(pROC))

script_arg <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if(length(script_arg) == 1) dirname(normalizePath(sub('^--file=', '', script_arg))) else '.'
source(file.path(script_dir, 'prediction_metrics.R'))

logf <- paste0(opt$out, '.log')
log_cat <- function(...) cat(..., file = logf, append = TRUE, sep = '')
cat('External validation (outer-fold model) started at ', as.character(start.time), '\nOptions:\n', file = logf, sep = '')
capture.output(print(opt), file = logf, append = TRUE)

###########
# Read the target cohort (same naming as Model_builder_V2_nested_wt_refit.R)
###########
predictors_list <- data.frame(fread(opt$predictors))
names(predictors_list)[1] <- 'predictor'
if(!is.null(predictors_list$group)) predictors_list$group <- gsub("[[:punct:]]", ".", predictors_list$group)
use_group_names <- opt$compare_predictors || (!is.null(predictors_list$group) && length(unique(predictors_list$group)) > 1)

Outcome <- fread(opt$pheno)
Outcome <- Outcome[complete.cases(Outcome),]
names(Outcome)[3] <- 'Outcome_var'
Outcome$IID <- paste0(Outcome$FID, ':', Outcome$IID)
Outcome$FID <- NULL
if(!is.na(opt$keep)){
	keep_file <- fread(opt$keep)
	names(keep_file)[1:2] <- c('FID','IID')
	Outcome <- Outcome[Outcome$IID %in% paste0(keep_file$FID, ':', keep_file$IID),]
}
binary_outcome <- length(unique(Outcome$Outcome_var)) == 2
log_cat('Target phenotype: ', nrow(Outcome), ' individuals; outcome is ', ifelse(binary_outcome, 'binary', 'continuous'), '.\n')

read_pred_file <- function(k){
	x <- fread(predictors_list$predictor[k])
	if(names(x)[1] == 'FID' & names(x)[2] == 'IID'){
		x$IID <- paste0(x$FID, ':', x$IID); x$FID <- NULL
	} else {
		names(x)[1] <- 'IID'; x$IID <- paste0(x$IID, ':', x$IID)
	}
	if(nrow(predictors_list) > 1){
		names(x)[-1] <- paste0('PredFile', k, '.', names(x)[-1])
		if(use_group_names) names(x)[-1] <- paste0('Group_', predictors_list$group[k], '.', names(x)[-1])
	} else if(opt$compare_predictors){
		names(x)[-1] <- paste0('Group_', gsub("[[:punct:]]", ".", names(x)[-1]), '.', names(x)[-1])
	}
	x
}
Predictors <- Reduce(function(a, b) merge(a, b, by = 'IID'), lapply(seq_len(nrow(predictors_list)), read_pred_file))
dat <- data.frame(merge(Outcome, Predictors, by = 'IID'))
names(dat) <- make.names(names(dat))   # data.frame() applies the same conversion in the development run
rm(Outcome, Predictors)

###########
# Load the outer-fold model
###########
f_new <- paste0(opt$model_prefix, '.', opt$group, '.fold', opt$fold, '.final_model.rds')   # Model_builder_V2_nested_wt_refit.R
f_old <- paste0(opt$model_prefix, opt$fold, '.final_model.rds')                              # Model_builder_V2_nested_wt.R
if(file.exists(f_new)){
	model_file <- f_new
} else if(file.exists(f_old)){
	model_file <- f_old
	log_cat('Note: using ', f_old, ' (earlier naming). In Model_builder_V2_nested_wt.R this file is overwritten by each predictor group, so it holds the last group fitted (All).\n')
} else {
	stop('Outer-fold model not found: ', f_new, ' or ', f_old)
}
model_label <- paste0('Fold', opt$fold)
models <- setNames(list(readRDS(model_file)), model_label)
log_cat('Model: ', model_label, ' = ', model_file, '\n')

# Predictor names expected by a model
model_vars <- function(fit){
	if(inherits(fit, 'glmnet')) rownames(fit$beta) else setdiff(names(coef(fit)), '(Intercept)')
}
needed <- unique(unlist(lapply(models, model_vars)))
needed <- needed[!grepl('^`?X?0`?$', needed)]   # constant column added for single-predictor glm models
missing <- setdiff(needed, names(dat))
if(length(missing) > 0) stop('Predictors used by the model are missing in the target cohort:\n  ', paste(missing, collapse = '\n  '),
	'\nCheck that the predictors file lists the same files, in the same order, with the same column names.')

# Complete cases on the predictors the model needs
n_before <- nrow(dat)
dat <- dat[complete.cases(dat[, c('Outcome_var', needed)]), ]
log_cat('Target cohort: ', nrow(dat), ' individuals with outcome and all model predictors (', n_before - nrow(dat), ' removed).\n')

predict_model <- function(fit){
	if(inherits(fit, 'glmnet')){
		x <- data.matrix(dat[, rownames(fit$beta), drop = FALSE])
		as.numeric(predict(fit, newx = x, s = fit$lambdaOpt, type = 'response'))
	} else {
		v <- setdiff(model_vars(fit), c('`0`', 'X0', '0'))
		tmp <- data.frame(0, dat[, v, drop = FALSE]); names(tmp)[1] <- '0'
		as.numeric(predict(fit, newdata = tmp, type = 'response'))
	}
}
preds <- lapply(models, predict_model)

###########
# Evaluate
###########
y <- dat$Outcome_var
model_is_binary <- inherits(models[[1]], 'lognet') || (inherits(models[[1]], 'glm') && family(models[[1]])$family == 'binomial')
label <- function(m) data.frame(Train = opt$train_label, Test = opt$test_label, Group = opt$group, Model = m)

eval_all <- NULL; thresh_all <- NULL; exp_all <- NULL; dec_all <- NULL; roc_list <- list(); dec_plot <- list()
for(m in names(preds)){
	p <- preds[[m]]
	if(binary_outcome){
		eval_all <- rbind(eval_all, cbind(label(m), binary_metrics(y, p, K = opt$outcome_pop_prev),
			if(model_is_binary) calibration_metrics(y, p) else NULL))
		roc_list[[m]] <- roc(as01(y), p, levels = c(0, 1), direction = '<', quiet = TRUE)
		if(model_is_binary){
			th <- threshold_metrics(y, p)
			thresh_all <- rbind(thresh_all, cbind(label(m), th))
			at05 <- th[th$Cutoff == 0.5, ]
			exp_all <- rbind(exp_all, cbind(label(m), Se = at05$Sensitivity, Sp = at05$Specificity,
				expected_ppv_npv(at05$Sensitivity, at05$Specificity, as.numeric(strsplit(opt$ppv_prev, ',')[[1]]))))
		}
	} else {
		# Continuous outcome: model trained for this outcome (external validation, Table S20),
		# or a diagnosis model applied to a severity outcome (Table S15; calibration columns not interpretable)
		eval_all <- rbind(eval_all, cbind(label(m), continuous_metrics(y, p)))
	}
	dec <- decile_table(y, p)
	dec_all <- rbind(dec_all, cbind(label(m), dec))
	dec_plot[[paste0('Elastic net (', m, ')')]] <- dec
}

write.table(eval_all, paste0(opt$out, '.ext_eval.txt'), sep = '\t', row.names = F, quote = F)
write.table(dec_all, paste0(opt$out, '.ext_deciles.txt'), sep = '\t', row.names = F, quote = F)
if(!is.null(thresh_all)) write.table(thresh_all, paste0(opt$out, '.ext_thresholds.txt'), sep = '\t', row.names = F, quote = F)
if(!is.null(exp_all)) write.table(exp_all, paste0(opt$out, '.ext_expected_ppv_npv.txt'), sep = '\t', row.names = F, quote = F)
if(length(roc_list) > 0) saveRDS(roc_list, paste0(opt$out, '.ext_roc.rds'))
pred_out <- data.frame(IID = dat$IID, Outcome = y, do.call(cbind, preds))
fwrite(pred_out, paste0(opt$out, '.ext_predictions.txt'), sep = '\t')
if(binary_outcome == model_is_binary){
	plot_calibration(dec_plot, paste0(opt$out, '.ext_calibration.pdf'),
		title = paste0('Trained in ', opt$train_label, ', tested in ', opt$test_label), binary = binary_outcome)
}

end.time <- Sys.time()
log_cat('Results saved with prefix ', opt$out, '.\nAnalysis finished at ', as.character(end.time),
	' (', round(as.numeric(difftime(end.time, start.time, units = 'mins')), 2), ' min).\n')
