# Cox proportional hazards model for incident mdd in UKB 

library(survival)

mdd_inc<-readRDS("mdd_onset_cox.rds")

cox_fit <- coxph(
  Surv(MHQ_years, MDD_onset) ~ .,
  data = mdd_inc[, c(
    "MDD_onset", "MHQ_years",
    "fh.depression", "fh.anxiety", "fh.mania", "fh.adhd",
    "fh.autism", "fh.eating_disorder", "fh.personality",
    "fh.psychosis", "fh.scz",
    "prs.MDDPGC3noUK", "prs.ANXnoUKB", "prs.ADHD", "prs.BIP",
    "prs.BMI", "prs.ASD", "prs.PTSD",
    "ChT", "sex",
    "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"
  )]
)

summary(cox_fit)

ph_test <- cox.zph(cox_fit)
print(ph_test)
plot(ph_test)

q()
