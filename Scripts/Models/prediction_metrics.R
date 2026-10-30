# prediction_metrics.R
# Shared evaluation functions for Model_builder_V2_nested_wt_refit.R (internal validation)
# and External_validation.R (external validation).
# Rujia Wang, 2026
#
# Metrics, matching the manuscript tables:
#   AUC with Hanley-McNeil 95% CI and DeLong 95% CI         (Tables S8, S20)
#   Observed and liability-scale R2 with 95% CI              (Tables S8, S15, S20)
#   Calibration: Brier score, calibration slope, intercept   (Table S14, Figure 3E-F)
#   Observed vs predicted by decile of predicted risk        (Figure 3E-F, Figure S5)
#   PPV, NPV, sensitivity, specificity at cut-offs 0.3-0.7   (Tables S16, S18)
#   Expected PPV/NPV at other prevalences                    (Table S17)

suppressMessages(library(pROC))

# Observed R2 -> liability-scale R2 (Lee et al., Genet Epidemiol 2012)
#   k: population prevalence; r2: observed-scale R2; p: proportion of cases in the sample
h2l_R2 <- function(k, r2, p) {
	x <- qnorm(1 - k)
	z <- dnorm(x)
	i <- z / k
	C <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
	theta <- i * ((p - k) / (1 - k)) * (i * ((p - k) / (1 - k)) - x)
	C * r2 / (1 + C * theta * r2)
}

# 95% CI for R2 using the large-sample approximation SE(R2) = sqrt(4 R2 (1 - R2)^2 / n)
# (n = number of cases for binary outcomes, total N for continuous outcomes, as in Table S20)
r2_ci <- function(r2, n) {
	se <- sqrt(4 * r2 * (1 - r2)^2 / n)
	c(low = max(0, r2 - 1.96 * se), high = min(1, r2 + 1.96 * se))
}

# 95% CI for AUC (Hanley & McNeil, Radiology 1982)
auc_ci_hanley <- function(auc, n_case, n_control) {
	q1 <- auc / (2 - auc)
	q2 <- 2 * auc^2 / (1 + auc)
	se <- sqrt((auc * (1 - auc) + (n_case - 1) * (q1 - auc^2) + (n_control - 1) * (q2 - auc^2)) / (n_case * n_control))
	c(low = auc - 1.96 * se, high = auc + 1.96 * se)
}

# Convert an outcome to 0/1 (accepts CONTROL/CASE factors, 0/1 or 1/2 codes)
as01 <- function(y) {
	if (is.factor(y) || is.character(y)) {
		y <- as.character(y)
		if (all(y %in% c('CONTROL', 'CASE'))) return(as.integer(y == 'CASE'))
		y <- as.numeric(y)
	}
	y <- as.numeric(y)
	if (all(y %in% c(1, 2))) y <- y - 1
	stopifnot(all(y %in% c(0, 1)))
	as.integer(y)
}

# Discrimination and explained variance for a binary outcome.
#   K: population prevalence for the liability-scale R2 (NA to skip)
binary_metrics <- function(y, p, K = NA) {
	y <- as01(y)
	n_case <- sum(y == 1); n_control <- sum(y == 0)
	roc_obj <- roc(y, p, levels = c(0, 1), direction = '<', quiet = TRUE)
	auc <- as.numeric(auc(roc_obj))
	ci_hm <- auc_ci_hanley(auc, n_case, n_control)
	ci_dl <- as.numeric(ci.auc(roc_obj, method = 'delong'))[c(1, 3)]
	# Standardised association (as in Model_builder: lm of scaled outcome on scaled score)
	lin <- summary(lm(scale(y) ~ scale(p)))
	R <- coef(lin)[2, 1]
	logit <- glm(y ~ scale(p), family = binomial)
	r2_obs <- R^2
	prop_case <- n_case / length(y)
	r2_liab <- if (!is.na(K)) h2l_R2(K, r2_obs, prop_case) else NA
	ci_obs <- r2_ci(r2_obs, n_case)
	# As in Table S20, the approximation is applied to the liability-scale R2 directly (n = number of cases)
	ci_liab <- if (!is.na(K)) r2_ci(r2_liab, n_case) else c(NA, NA)
	data.frame(
		N = length(y), Ncase = n_case, Ncont = n_control, Case_prop = prop_case,
		R = R, R_SE = coef(lin)[2, 2], P = coef(lin)[2, 4],
		OR_perSD = exp(coef(logit)[2]),
		AUC = auc, AUC_low_HM = ci_hm[1], AUC_high_HM = ci_hm[2],
		AUC_low_DeLong = ci_dl[1], AUC_high_DeLong = ci_dl[2],
		R2_obs = r2_obs, R2_obs_low = ci_obs[1], R2_obs_high = ci_obs[2],
		K = K, R2_liab = r2_liab, R2_liab_low = ci_liab[1], R2_liab_high = ci_liab[2],
		row.names = NULL)
}

# Calibration of predicted probabilities against a binary outcome.
#   cal_intercept / cal_slope: logistic recalibration model y ~ logit(p), as used for XGBoost (xgb_mdd.R)
#     and reported in Table S14 (ideal: intercept 0, slope 1)
#   cal_in_large: calibration-in-the-large, intercept of y ~ offset(logit(p)) (ideal 0;
#     negative = risk overestimated, positive = risk underestimated)
#   brier: mean squared difference between predicted probability and outcome
calibration_metrics <- function(y, p) {
	y <- as01(y)
	eps <- 1e-6
	p <- pmin(pmax(p, eps), 1 - eps)
	lp <- qlogis(p)
	fit <- glm(y ~ lp, family = binomial)
	fit_large <- glm(y ~ offset(lp), family = binomial)
	data.frame(
		Brier = mean((p - y)^2),
		Cal_slope = unname(coef(fit)[2]),
		Cal_intercept = unname(coef(fit)[1]),
		Cal_in_large = unname(coef(fit_large)[1]),
		Mean_pred = mean(p), Obs_prop = mean(y),
		row.names = NULL)
}

# Observed outcome against predicted score by quantile group (deciles by default).
#   Binary outcome: observed proportion with 95% Wilson CI vs mean predicted probability (calibration plot).
#   Continuous outcome: observed mean with SD and 95% CI vs mean predicted score (Figure S5).
decile_table <- function(y, p, n_groups = 10) {
	binary <- length(unique(y)) == 2
	if (binary) y <- as01(y) else y <- as.numeric(y)
	brks <- unique(quantile(p, probs = seq(0, 1, length.out = n_groups + 1), na.rm = TRUE, type = 7))
	grp <- cut(p, breaks = brks, include.lowest = TRUE, labels = FALSE)
	out <- do.call(rbind, lapply(sort(unique(grp)), function(g) {
		yy <- y[grp == g]; pp <- p[grp == g]; n <- length(yy)
		if (binary) {
			ph <- mean(yy); z <- 1.96
			centre <- (ph + z^2 / (2 * n)) / (1 + z^2 / n)
			half <- z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
			data.frame(Group = g, N = n, Mean_pred = mean(pp), Observed = ph,
				Obs_SD = NA, Obs_low = centre - half, Obs_high = centre + half)
		} else {
			s <- sd(yy)
			data.frame(Group = g, N = n, Mean_pred = mean(pp), Observed = mean(yy),
				Obs_SD = s, Obs_low = mean(yy) - 1.96 * s / sqrt(n), Obs_high = mean(yy) + 1.96 * s / sqrt(n))
		}
	}))
	out
}

# PPV, NPV, sensitivity and specificity at a set of probability cut-offs (Tables S16, S18)
threshold_metrics <- function(y, p, cutoffs = c(0.3, 0.4, 0.5, 0.6, 0.7)) {
	y <- as01(y)
	do.call(rbind, lapply(cutoffs, function(cut) {
		pos <- p >= cut
		TP <- sum(pos & y == 1); FP <- sum(pos & y == 0)
		TN <- sum(!pos & y == 0); FN <- sum(!pos & y == 1)
		data.frame(Cutoff = cut, TP = TP, FP = FP, TN = TN, FN = FN,
			PPV = TP / (TP + FP), NPV = TN / (TN + FN),
			Sensitivity = TP / (TP + FN), Specificity = TN / (TN + FP))
	}))
}

# Expected PPV and NPV at other prevalences from observed sensitivity and specificity (Table S17)
expected_ppv_npv <- function(se, sp, prev = c(0.10, 0.20, 0.36, 0.69)) {
	data.frame(Prevalence = prev,
		PPV = se * prev / (se * prev + (1 - sp) * (1 - prev)),
		NPV = sp * (1 - prev) / (sp * (1 - prev) + (1 - se) * prev))
}

# Explained variance and calibration-in-the-mean for a continuous outcome
#   (number of episodes, age of onset; Tables S15 and S20).
#   R / R2_obs: standardised association between score and outcome.
#   Cal_slope / Cal_intercept: lm(outcome ~ score); meaningful only when the score is on the outcome scale
#   (a model trained for that outcome); not interpretable when an MDD-diagnosis model is applied (Table S15).
continuous_metrics <- function(y, p) {
	y <- as.numeric(y)
	lin <- summary(lm(scale(y) ~ scale(p)))
	R <- coef(lin)[2, 1]
	ci <- r2_ci(R^2, length(y))
	raw <- coef(lm(y ~ p))
	data.frame(N = length(y), R = R, R_SE = coef(lin)[2, 2], P = coef(lin)[2, 4],
		R2_obs = R^2, R2_obs_low = ci[1], R2_obs_high = ci[2],
		Cal_slope = unname(raw[2]), Cal_intercept = unname(raw[1]),
		MSE = mean((y - p)^2), Mean_pred = mean(p), Mean_obs = mean(y),
		row.names = NULL)
}

# Calibration plot from one or more decile tables (list names are used as legend labels)
plot_calibration <- function(dec_list, file, title = '', binary = TRUE) {
	suppressMessages(library(ggplot2))
	d <- do.call(rbind, lapply(names(dec_list), function(nm) cbind(Model = nm, dec_list[[nm]])))
	lims <- range(c(d$Mean_pred, d$Obs_low, d$Obs_high), na.rm = TRUE)
	if (binary) lims <- c(0, 1)
	g <- ggplot(d, aes(x = Mean_pred, y = Observed, colour = Model)) +
		geom_abline(intercept = 0, slope = 1, linetype = 'dashed', colour = 'grey50') +
		geom_errorbar(aes(ymin = Obs_low, ymax = Obs_high), width = 0) +
		geom_line() + geom_point() +
		coord_equal(xlim = lims, ylim = lims) +
		labs(x = if (binary) 'Mean predicted probability' else 'Mean predicted score',
			y = if (binary) 'Observed proportion of MDD (95% CI)' else 'Observed mean (95% CI)',
			title = title, colour = NULL) +
		theme_minimal(base_size = 12) + theme(legend.position = 'bottom')
	dev_fun <- if (capabilities('cairo')) cairo_pdf else pdf
	ggsave(file, g, width = 6, height = 6.4, device = dev_fun, bg = 'white')
	invisible(g)
}
