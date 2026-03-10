#### MDD episodes in UKB
#### Rujia Wang
#### 2024-01-08

library(ukbkings)
library(dplyr)
library(stringr)
library(data.table)

project_dir <- "/datasets/ukbiobank/ukb82087"

#--------------------------------------------------
# 1. Extract UKB fields
#--------------------------------------------------

f <- bio_field(project_dir)

f %>%
  select(field, name) %>%
  filter(str_detect(field, "29033|29011|29012|20442|20441|20446|4620")) %>%
  bio_field_add("MDD_episodes.txt")

bio_phen(
  project_dir,
  field = "MDD_episodes.txt",
  out = "MDD_episodes"
)

mydata <- readRDS("MDD_episodes.rds")

#--------------------------------------------------
# 2. Helper function to recode episode counts
#--------------------------------------------------

cap_episodes <- function(x) {
  case_when(
    x > 0 & x < 13 ~ as.numeric(x),
    x > 12 ~ 13,
    TRUE ~ NA_real_
  )
}

#--------------------------------------------------
# 3. Recode MDD episode counts from each visit
#--------------------------------------------------
# 4620-0.0 = baseline (2006-2010)
# 4620-1.0 = first repeat assessment (2012-2013)
# 4620-2.0 = imaging visit (2014+)
# 4620-3.0 = repeat imaging visit (2019+)
# 20442-0.0, 29033-0.0 = MHQ-derived episode fields

mydata <- mydata %>%
  mutate(
    MDD_episodes_v1   = cap_episodes(`4620-0.0`),
    MDD_episodes_v2   = cap_episodes(`4620-1.0`),
    MDD_episodes_v3   = cap_episodes(`4620-2.0`),
    MDD_episodes_v4   = cap_episodes(`4620-3.0`),
    MDD_episodes_MHQ1 = cap_episodes(`20442-0.0`),
    MDD_episodes_MHQ2 = cap_episodes(`29033-0.0`)
  )

#--------------------------------------------------
# 4. Combine MHQ episode definitions
#--------------------------------------------------
# Take the larger reported count across MHQ1 and MHQ2

mydata <- mydata %>%
  mutate(
    MDD_episodes = pmax(MDD_episodes_MHQ1, MDD_episodes_MHQ2, na.rm = TRUE),
    MDD_episodes = if_else(is.infinite(MDD_episodes), NA_real_, MDD_episodes)
  )

#--------------------------------------------------
# 5. Combine MHQ and repeated-measure episode counts
#--------------------------------------------------
# Sequentially retain the maximum across MHQ + visits 1-4

mydata <- mydata %>%
  mutate(
    MDD_episodes_MHQ_v1   = pmax(MDD_episodes, MDD_episodes_v1, na.rm = TRUE),
    MDD_episodes_MHQ_v1   = if_else(is.infinite(MDD_episodes_MHQ_v1), NA_real_, MDD_episodes_MHQ_v1),
    
    MDD_episodes_MHQ_v12  = pmax(MDD_episodes_MHQ_v1, MDD_episodes_v2, na.rm = TRUE),
    MDD_episodes_MHQ_v12  = if_else(is.infinite(MDD_episodes_MHQ_v12), NA_real_, MDD_episodes_MHQ_v12),
    
    MDD_episodes_MHQ_v123 = pmax(MDD_episodes_MHQ_v12, MDD_episodes_v3, na.rm = TRUE),
    MDD_episodes_MHQ_v123 = if_else(is.infinite(MDD_episodes_MHQ_v123), NA_real_, MDD_episodes_MHQ_v123),
    
    MDD_episodes_MHQ_v1234 = pmax(MDD_episodes_MHQ_v123, MDD_episodes_v4, na.rm = TRUE),
    MDD_episodes_MHQ_v1234 = if_else(is.infinite(MDD_episodes_MHQ_v1234), NA_real_, MDD_episodes_MHQ_v1234)
  )

#--------------------------------------------------
# 6. Quick QC checks
#--------------------------------------------------

table(mydata$MDD_episodes, mydata$MDD, useNA = "ifany")
table(mydata$MDD_episodes_MHQ_v1234, useNA = "ifany")

#--------------------------------------------------
# 7. Merge with broad MDD phenotype
#--------------------------------------------------

mdd_broad <- fread("mdd_broad_control_rmANX_covariates_340062.txt")

tep <- mdd_broad %>%
  filter(sample == "UKB") %>%
  select(FID, MDD_broad)

mydata <- mydata %>%
  left_join(tep, by = c("eid" = "FID"))

#--------------------------------------------------
# 8. Define GWAS phenotypes
#--------------------------------------------------
# For GWAS:
# - cases: keep episode count
# - controls: set to 0
# - controls with non-missing repeated-episode information: set to NA

mydata <- mydata %>%
  mutate(
    MDD_episodes_mhq = case_when(
      MDD_broad == 0 & !is.na(MDD_episodes_MHQ_v1234) ~ NA_real_,
      MDD_broad == 1 ~ MDD_episodes,
      MDD_broad == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    MDD_episodes_mhq_4R = case_when(
      MDD_broad == 0 & !is.na(MDD_episodes_MHQ_v1234) ~ NA_real_,
      MDD_broad == 1 ~ MDD_episodes_MHQ_v1234,
      MDD_broad == 0 ~ 0,
      TRUE ~ NA_real_
    )
  )

# optional cleaned diagnosis variable
mydata <- mydata %>%
  mutate(
    MDD_con_rm_episodes = case_when(
      MDD_broad == 0 & !is.na(MDD_episodes_MHQ_v1234) ~ NA_real_,
      TRUE ~ as.numeric(MDD)
    )
  )

#--------------------------------------------------
# 9. Save output
#--------------------------------------------------

save(mydata, file = "MDD_episodes.Rdata")

q()
