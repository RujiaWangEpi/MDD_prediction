########################################################################
#### Figure4. Ridgeline plot of prediction score among MDD episodes ####
#### Author   Rujia Wang                                            ####
#### Email    rujia.1.wang@kcl.ac.uk                                ####
#### Date     2025-02-25                                            #### 
########################################################################

library(tidyverse)
library(ggplot2)
library(readxl)
library(ggridges)
library(dplyr)
library(ggrepel)
library(patchwork)

setwd("working path")

tep<-readRDS("UKB_mdd_pred_across_episodes.rds")
tep1<-readRDS("GLAD_mdd_pred_across_episodes.rds")

tep1<-subset(tep1,!is.na(tep1$mdd_episode),select = c("obs","pred","mdd_episode"))
tep1$Cohort<-"GLAD"

tep$Cohort<-"UKB"
tep<-rbind(tep1,tep3)

order_levels <- c("Control","Single_episode", "2_episodes", "3_episodes","4_episodes","5_episodes", "6_episodes","713_episodes")
order_labels <- c("Control","Single episode", "2 episodes", "3 episodes","4 episodes","5 episodes", "6 episodes", "≥7 episodes")

tep$mdd_episode <- factor(tep$mdd_episode, levels = order_levels, labels = order_labels)

cohort_order_levels <- c("GLAD", "UKB")
cohort_order_labels <- c("GLAD+", "UKB")

tep$Cohort<-factor(tep$Cohort, levels = cohort_order_levels, labels = cohort_order_labels)
#tep$Cohort<-factor(tep$Cohort, levels = rev(cohort_order_levels), labels = rev(cohort_order_labels))

plot <- tep %>%
  ggplot() +
  aes(x = pred, 
      y = mdd_episode,
      fill = mdd_episode) +
  geom_density_ridges2(scale = 5, bandwidth = 0.0671, height = 0.1,
                       jittered_points = TRUE, point_alpha = 0.001) + 
  scale_fill_manual(values = c("#DEF1D0","#AFE5C4","#62C5BC","#28A4B3","#068FAB","#0078A2","#1C518F","#2D3184")) +
  scale_color_viridis_d(alpha = 0.9) +  
  labs(x = "Liability of being MDD cases",
       y = "MDD episodes",
       fill = "MDD episodes") +  
  theme_minimal() +
  theme(axis.title = element_text(size = 16),  
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 14),  
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_cartesian(xlim = c(0, 1)) +
  facet_wrap(~ Cohort, scales = "free_y", ncol = 2)

ggsave("working path/Figures/Figure4.UKB_GLAD_MDD_episodes_OR_ridgeline_20250222.png", plot, width = 16, height = 7, units = "in", dpi = 600)

q()

