# MDD_prediction
"**Combining family history, polygenic scores, and childhood trauma to predict major depressive disorder**"

Wang R, Davies H, Lee SY, Zvrskovec J, Hübel C, et al. Preprint: https://www.medrxiv.org/content/10.1101/2025.09.22.25336356v1

Major depressive disorder (MDD) is one of the most common mental health conditions worldwide. However, predicting who will develop depression — or who is at higher risk — remains difficult. This is partly because depression arises from a complex combination of genetic, family, environmental, and social factors.

In this project, we examine whether combining several types of risk information can improve prediction of depression. These include:

1) Family history of psychiatric disorders

2) Genetic risk measured using polygenic risk scores (PRS)

3) Exposure to childhood trauma

4) Demographic factors such as sex and genetic principal components (PCs)

We analyzed data from two large cohorts: the GLAD+ study, which focuses on individuals with anxiety and depression, and the UK Biobank, a large population-based study.

Using machine learning and statistical models, we evaluated how well these factors predict:

1) Whether someone has depression

2) The number of depressive episodes and age of onset (severity)

Our results show that combining multiple risk factors improves prediction compared with using genetic information alone. In particular, family history and childhood trauma were consistently stronger predictors of depression than genetic risk scores. While genetic information contributes useful insights, traditional clinical risk factors remain highly informative.

These findings highlight that integrating familial, environmental, and genetic information provides a more complete picture of depression risk. The work contributes to ongoing efforts to develop more accurate and clinically meaningful approaches for identifying individuals at elevated risk of depression.

## Repository contents

| Folder | Contents |
|---|---|
| `Scripts/Pheno/` | UK Biobank phenotype definitions (number of episodes, age of onset, family history from MHQ2, childhood trauma). Lifetime MDD in MHQ1/MHQ2: https://github.com/ColemanResearchGroup/MHQ2. GLAD+ phenotypes were derived with the same instruments from the COPING survey (Supplementary Methods). |
| `Scripts/Models/` | Nested elastic net models, full-sample refits, external validation, evaluation metrics (AUC, liability R², calibration, PPV/NPV), DeLong comparisons and the XGBoost sensitivity analysis. `Scripts/Models/README.md` maps each script to the manuscript tables and figures. |

Polygenic scores were computed with the GenoPred pipeline (https://github.com/opain/GenoPred; MegaPRS); settings are given in the Supplementary Methods.
