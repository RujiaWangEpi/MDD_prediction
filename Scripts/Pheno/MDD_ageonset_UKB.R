#### MDD age at onset in UKB ####
#### Rujia Wang ####
#### 2024-04-23 ####

library(ukbkings)
library(dplyr)
library(stringr)

project_dir <- "/datasets/ukbiobank/ukb82087"

#--------------------------------------------------
# 1. Extract UKB fields
#--------------------------------------------------

f <- bio_field(project_dir)

head(f)
glimpse(f)

f %>%
  distinct(basket)

f %>%
  select(field, name) %>%
  filter(str_detect(field, "29034|20433")) %>%
  bio_field_add("MDD_ageonset.txt")

bio_phen(
  project_dir,
  field = "MDD_ageonset.txt",
  out = "MDD_ageonset"
)

MDD_ageonset <- readRDS("MDD_ageonset.rds")

#--------------------------------------------------
# 2. Helper function to clean age-at-onset values
#--------------------------------------------------
# Negative values are UKB missing / invalid codes and are set to NA

clean_ageonset <- function(x) {
  case_when(
    x < 0 ~ NA_real_,
    TRUE ~ as.numeric(x)
  )
}

#--------------------------------------------------
# 3. Derive MHQ1 and MHQ2 age-at-onset variables
#--------------------------------------------------

MDD_ageonset <- MDD_ageonset %>%
  mutate(
    mdd_ageonset_MHQ1 = clean_ageonset(`20433-0.0`),
    mdd_ageonset_MHQ2 = clean_ageonset(`29034-0.0`)
  )

#--------------------------------------------------
# 4. Define final age at onset
#--------------------------------------------------
# Use the earliest non-missing age at onset reported across MHQ1 and MHQ2

MDD_ageonset <- MDD_ageonset %>%
  mutate(
    mdd_ageonset = pmin(mdd_ageonset_MHQ1, mdd_ageonset_MHQ2, na.rm = TRUE),
    mdd_ageonset = if_else(is.infinite(mdd_ageonset), NA_real_, mdd_ageonset)
  )

#--------------------------------------------------
# 5. QC checks
#--------------------------------------------------

table(MDD_ageonset$mdd_ageonset_MHQ1, useNA = "ifany")
table(MDD_ageonset$mdd_ageonset_MHQ2, useNA = "ifany")
table(MDD_ageonset$mdd_ageonset, useNA = "ifany")

#--------------------------------------------------
# 6. Save output
#--------------------------------------------------

save(MDD_ageonset, file = "MDD_ageonset.Rdata")

q()
