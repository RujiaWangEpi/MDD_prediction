#######################################################################
#### Figure3.The effect of group and individual predictors for MDD ####
#### Author   Rujia Wang                                           ####
#### Email    rujia.1.wang@kcl.ac.uk                               ####
#### Date     2025-02-25                                           #### 
#######################################################################

library(tidyverse)
library(ggplot2)
library(data.table)
library(readxl)
library(ggridges)
library(dplyr)
library(ggrepel)
library(patchwork)

setwd("working path")

df <- tribble(
  ~predictors, ~liability, ~cohort,
  "mFH", 16.63, "GLAD+",
  "mPRS", 7.33, "GLAD+",
  "ChT", 10.64, "GLAD+",
  "Demo", 9.77, "GLAD+",
  "All", 32.91, "GLAD+",
  "mFH", 12.66, "UKB",
  "mPRS", 3.68, "UKB",
  "ChT", 7.34, "UKB",
  "Demo", 6.05, "UKB",
  "All", 23.08, "UKB"
)

predictor_order <- c("mFH", "mPRS", "ChT", "Demo", "All")

# Convert 'predictors' to factor with specified levels/order
df$predictors <- factor(df$predictors, levels = predictor_order)

# Plotting using ggplot2

predictor_colors <- c("#0099ff", "#4daf4a", "#FF5A5F", "#ffa500", "#7664AE")

##### Figure 3: R2 and AUC #########
# Create a PNG device for saving the combined plot
png("/output/GLAD_UKB_MDD_combined_plot_r2_auc.png", width = 2400, height = 2400, res = 300)

# Set up the layout for the plot
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2) + 0.1, oma = c(0, 0, 2, 0)) 

# Predictor colors
predictor_colors <- c("#0099ff", "#4daf4a", "#FF5A5F", "#ffa500", "#7664AE")

# Read in ROC data for GLAD
roc_model1 <- readRDS("GLAD/mdd_all_pred.Indep_auc_fh.rds")
roc_model2 <- readRDS("GLAD/mdd_all_pred.Indep_auc_prs.rds")
roc_model3 <- readRDS("GLAD/mdd_all_pred.Indep_auc_ct.rds")
roc_model4 <- readRDS("GLAD/mdd_all_pred.Indep_auc_demo.rds")
roc_model5 <- readRDS("GLAD/mdd_all_pred.Indep_auc_all.rds")

# Read in ROC data for UKB
roc_model6 <- readRDS("UKB/nested_ukb_mdd_allpred.Indep_auc_fh.rds")
roc_model7 <- readRDS("UKB/nested_ukb_mdd_allpred.Indep_auc_prs.rds")
roc_model8 <- readRDS("UKB/nested_ukb_mdd_allpred.Indep_auc_ct.rds")
roc_model9 <- readRDS("UKB/nested_ukb_mdd_allpred.Indep_auc_demo.rds")
roc_model10 <- readRDS("UKB/nested_ukb_mdd_allpred.Indep_auc_all.rds")

### First plot: GLAD Liability ###
glad_data <- subset(df, cohort == "GLAD+")
barplot(glad_data$liability, names.arg = glad_data$predictors, col = predictor_colors, ylim = c(0, 36), 
        main = "A. GLAD+ Liability_R2", ylab = "Liability_R2 (%)")
text(x = barplot(glad_data$liability, plot = FALSE), y = glad_data$liability, label = sprintf("%.2f", glad_data$liability), pos = 3, cex = 1.0, col = "black")

### Second plot: UKB Liability ###
ukb_data <- subset(df, cohort == "UKB")
barplot(ukb_data$liability, names.arg = ukb_data$predictors, col = predictor_colors, ylim = c(0, 36), 
        main = "B. UKB Liability_R2", ylab = "Liability_R2 (%)")
text(x = barplot(ukb_data$liability, plot = FALSE), y = ukb_data$liability, label = sprintf("%.2f", ukb_data$liability), pos = 3, cex = 1.0, col = "black")

### Third plot: GLAD ROC ###
plot(roc_model1, legacy.axes = TRUE, percent = TRUE, xlab = "False-positive rate", ylab = "True-positive rate", col = "#0099ff", lwd = 2, cex.lab = 1.2)
plot(roc_model2, col = "#4daf4a", add = TRUE)
plot(roc_model3, col = "#FF5A5F", add = TRUE)
plot(roc_model4, col = "#ffa500", add = TRUE)
plot(roc_model5, col = "#481B6D", add = TRUE)
title(main = "C. GLAD+ ROC", adj = 0.5, cex.main = 1.2, line = 1)  
legend("bottomright", legend = c("mFH_AUC: 0.74", "mPRS_AUC: 0.67", "ChT_AUC: 0.67", "Demo_AUC: 0.67", "All_AUC: 0.84"),
       col = c("#0099ff", "#4daf4a", "#FF5A5F", "#ffa500", "#481B6D"), lwd = 2, cex = 1.0)

### Fourth plot: UKB ROC ###
plot(roc_model6, legacy.axes = TRUE, percent = TRUE, xlab = "False-positive rate", ylab = "True-positive rate", col = "#0099ff", lwd = 2, cex.lab = 1.2)
plot(roc_model7, col = "#4daf4a", add = TRUE)
plot(roc_model8, col = "#FF5A5F", add = TRUE)
plot(roc_model9, col = "#ffa500", add = TRUE)
plot(roc_model10, col = "#481B6D", add = TRUE)
title(main = "D. UKB ROC", adj = 0.5, cex.main = 1.2, line = 1) 
legend("bottomright", legend = c("mFH_AUC: 0.65", "mPRS_AUC: 0.60", "ChT_AUC: 0.61", "Demo_AUC: 0.62", "All_AUC: 0.74"),
       col = c("#0099ff", "#4daf4a", "#FF5A5F", "#ffa500", "#481B6D"), lwd = 2, cex = 1.0)

# Close the PNG device
dev.off()

#### Figure 3E and 3F: Individual predictors ####

# Read in the GLAD data
# tep <- read_excel("GLAD_MDD_individual_predictors_20241012.xlsx")
# Read in the GLAD map UKB data 
tep <- read_excel("GLAD_MDD_individual_predictors_map_UKB_20241012.xlsx")
tep$Group <- as.factor(tep$Group)
tep$Z_score <- tep$R / tep$SE

# Read in the UKB data
tep1 <- read_excel("UKB_MDD_individual_predictors_20241120.xlsx")
tep1$Group <- as.factor(tep1$Group)
tep1$Z_score <- tep1$R / tep1$SE

# Open PNG device for E-F (1x2 layout)
png("Figure3_GLAD_UKB_MDD_Individual_Predictors.png", 
    width = 2400, height = 1200, res = 300)

# Create E: Individual Predictors in GLAD
plot1 <- ggplot(tep, aes(x = log2(OR), y = Z_score, color = Group, label = Model)) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = ifelse(log2(OR) > 0.75, Model, "")), size = 3.5, 
                  box.padding = unit(0.15, "lines"), nudge_x = 0.3, nudge_y = 0.2, segment.color = NA) +
  labs(x = "log2(OR)", y = "Z-score", title = "E. Individual predictors in GLAD+") +
  theme_minimal() +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("#F5546E","#FFB400", "#481B6D", "#00A04B")) +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
        legend.position = "none",  # Remove legend
        panel.border = element_rect(color = "black", fill = NA, size = 0.7),  
        panel.grid = element_blank(), # Remove grid lines
        plot.margin = margin(15, 20, 15, 10)) + xlim(-0.1,2.0)  

# Create F: Individual Predictors in UKB
plot2 <- ggplot(tep1, aes(x = log2(OR), y = Z_score, color = Group, label = Model)) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = ifelse(Z_score > 50, Model, "")), size = 3.5, 
                  box.padding = unit(0.15, "lines"), nudge_x = 0.3, nudge_y = 0.2, segment.color = NA) +
  labs(x = "log2(OR)", y = "Z-score", title = "F. Individual predictors in UKB") +
  theme_minimal() +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("#F5546E","#FFB400", "#481B6D", "#00A04B")) +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
        legend.position = "none",  # Remove legend
        panel.border = element_rect(color = "black", fill = NA, size = 0.7), 
        panel.grid = element_blank(),
        plot.margin = margin(15, 20, 15, 10)) + xlim(-0.05,1.5)  

# Combine the plots using grid.arrange()
grid.arrange(plot1, plot2, ncol = 2, widths = c(1, 1), padding = unit(10, "cm"))  

dev.off() 

q()

