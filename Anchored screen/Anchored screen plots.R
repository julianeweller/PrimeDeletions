library(tidyverse)
library(dplyr)
library(ggpointdensity)
library(viridis)
library(multcomp)
library(egg)
library(grid)
library(ineq)
library(scales)
library(reshape2)
library(janitor)

theme_pie <-   theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        text = element_text(family = "Helvetica"),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        plot.background = element_blank())

colours_dip <- c("#F9B233","#BE1622","#1D71B8","#A0D2F3","#006633","#7BA68F","#82368C","#AC8BBF","#F39200","#F9CDAB","#B2B2B2", "#4C5454","#E6007E", "#DEDC00") #Deletion # Homology Arm, # guides, 

path = "/"
outdir = paste0(path, "figures/")
outdir_supp = paste0(path, "supplements/")

# Load data
dipper <- read_csv(paste0(path, "data.csv"))


##########################################################################
################################ FIGURE 1 ################################
##########################################################################

####### Replicate correlation
interim_data <- read_csv(paste0(path, "files/Dipper_replicatecorrelation.csv")) %>% 
  mutate(L2FC = as.numeric(L2FC)) %>%
  filter(!is.na(L2FC))  %>% 
  dplyr::select(ID_replicate, unique, L2FC)

reshaped_data <- interim_data %>%
  pivot_wider(
    names_from = ID_replicate,
    values_from = L2FC
  )

reshaped_data <- reshaped_data %>%
  clean_names()

corr_replicates <- list(
  HEK_DPr04_1B_vs_HEK_DPr04_2B = cor(reshaped_data$x51_1_hek_d_pr04_1b, reshaped_data$x51_1_hek_d_pr04_2b, use = "complete.obs", method = "pearson"),
  HEK_DPr04_1C_vs_HEK_DPr04_2C = cor(reshaped_data$x51_1_hek_d_pr04_1c, reshaped_data$x51_1_hek_d_pr04_2c, use = "complete.obs", method = "pearson"),
  HEK_BGHr01_1_vs_HEK_BGHr01_2 = cor(reshaped_data$x51_2_hek_bg_hr01_1, reshaped_data$x51_2_hek_bg_hr01_2, use = "complete.obs", method = "pearson"),
  HAP_BGHr01_3_vs_HAP_BGHr01_4 = cor(reshaped_data$x51_2_hap_bg_hr01_3, reshaped_data$x51_2_hap_bg_hr01_4, use = "complete.obs", method = "pearson"),
  HEK_DPr04_1a_vs_HEK_DPr04_1b = cor(reshaped_data$x51_3_hek_d_pr04_1a, reshaped_data$x51_3_hek_d_pr04_2a, use = "complete.obs", method = "pearson"),
  HEK_DPr08_1a_vs_HEK_DPr08_1b = cor(reshaped_data$x51_3_hek_d_pr08_1b, reshaped_data$x51_3_hek_d_pr08_2b, use = "complete.obs", method = "pearson"),
  HEK_BGH_1a_vs_HEK_BGH_1b = cor(reshaped_data$x51_3_hek_bg_hr01_1c, reshaped_data$x51_3_hek_bg_hr01_2c, use = "complete.obs", method = "pearson"),
  HAP_DPr04_3a_vs_HAP_DPr04_4b = cor(reshaped_data$x51_3_hap_d_pr04_3a, reshaped_data$x51_3_hap_d_pr04_4a, use = "complete.obs", method = "pearson"),
  HAP_DPr08_3a_vs_HAP_DPr08_4b = cor(reshaped_data$x51_3_hap_d_pr08_3b, reshaped_data$x51_3_hap_d_pr08_4b, use = "complete.obs", method = "pearson"),
  HAP_BGH_3a_vs_HAP_BGH_4b = cor(reshaped_data$x51_3_hap_bg_hr01_3c, reshaped_data$x51_3_hap_bg_hr01_4c, use = "complete.obs", method = "pearson")
)

# Correlation between cell lines: L2FC_HAP_BGHr01 vs L2FC_HEK_BGHr01, L2FC_HAP_DPr04 vs L2FC_HEK_DPr04, L2FC_HAP_DPr08 vs L2FC_HEK_DPr08
corr_sameanchor_diffcellline <- list(
  HAP_vs_HEK_BGHr01 = cor(dipper_all$L2FC_HAP_BGHr01, dipper_all$L2FC_HEK_BGHr01, use = "complete.obs", method = "pearson"),
  HAP_vs_HEK_DPr04 = cor(dipper_all$L2FC_HAP_DPr04, dipper_all$L2FC_HEK_DPr04, use = "complete.obs", method = "pearson"),
  HAP_vs_HEK_DPr08 = cor(dipper_all$L2FC_HAP_DPr08, dipper_all$L2FC_HEK_DPr08, use = "complete.obs", method = "pearson"))

# Correlation between anchors: L2FC_HAP_BGHr01 vs L2FC_HAP_DPr04, L2FC_HAP_BGHr01 vs L2FC_HAP_DPr08, L2FC_HEK_BGHr01 vs L2FC_HEK_DPr04, L2FC_HEK_BGHr01 vs L2FC_HEK_DPr08
corr_diffanchors_samecellline <- list(
  HAP_BGHr01_vs_DPr04 = cor(dipper_all$L2FC_HAP_BGHr01, dipper_all$L2FC_HAP_DPr04, use = "complete.obs", method = "pearson"),
  HAP_BGHr01_vs_DPr08 = cor(dipper_all$L2FC_HAP_BGHr01, dipper_all$L2FC_HAP_DPr08, use = "complete.obs", method = "pearson"),
  HAP_DPr08_vs_DPr04 = cor(dipper_all$L2FC_HAP_DPr08, dipper_all$L2FC_HAP_DPr04, use = "complete.obs", method = "pearson"),
  HEK_BGHr01_vs_DPr04 = cor(dipper_all$L2FC_HEK_BGHr01, dipper_all$L2FC_HEK_DPr04, use = "complete.obs", method = "pearson"),
  HEK_BGHr01_vs_DPr08 = cor(dipper_all$L2FC_HEK_BGHr01, dipper_all$L2FC_HEK_DPr08, use = "complete.obs", method = "pearson"),
  HEK_DPr08_vs_DPr04 = cor(dipper_all$L2FC_HEK_DPr08, dipper_all$L2FC_HEK_DPr04, use = "complete.obs", method = "pearson"))

corr_diffanchors_diffcellline <- list(
  HAP_BGHr01_vs_HEK_DPr04 = cor(dipper_all$L2FC_HAP_BGHr01, dipper_all$L2FC_HEK_DPr04, use = "complete.obs", method = "pearson"),
  HAP_BGHr01_vs_HEK_DPr08 = cor(dipper_all$L2FC_HAP_BGHr01, dipper_all$L2FC_HEK_DPr08, use = "complete.obs", method = "pearson"),
  HAP_DPr08_vs_HEK_DPr04 = cor(dipper_all$L2FC_HAP_DPr08, dipper_all$L2FC_HEK_DPr04, use = "complete.obs", method = "pearson"),
  HAP_DPr08_vs_HEK_BGHr01 = cor(dipper_all$L2FC_HAP_DPr08, dipper_all$L2FC_HEK_BGHr01, use = "complete.obs", method = "pearson"),
  HAP_DPr04_vs_HEK_DPr08 = cor(dipper_all$L2FC_HAP_DPr04, dipper_all$L2FC_HEK_DPr08, use = "complete.obs", method = "pearson"),
  HAP_DPr04_vs_HEK_BGHr01 = cor(dipper_all$L2FC_HAP_DPr04, dipper_all$L2FC_HEK_BGHr01, use = "complete.obs", method = "pearson"))

corr_meananchors <- list(
  BGHr01_vs_DPr04 = cor(dipper_all$L2FC_BGHr01, dipper_all$L2FC_DPr04, use = "complete.obs", method = "pearson"),
  BGHr01_vs_DPr08 = cor(dipper_all$L2FC_BGHr01, dipper_all$L2FC_DPr08, use = "complete.obs", method = "pearson"),
  DPr08_vs_DPr04 = cor(dipper_all$L2FC_DPr08, dipper_all$L2FC_DPr04, use = "complete.obs", method = "pearson"))

corr_meancelllines <- list(
  HEK_vs_HAP = cor(dipper_all$L2FC_HEK, dipper_all$L2FC_HAP, use = "complete.obs", method = "pearson"))

df_replicates <- data.frame(List = "Replicates", Correlation = unlist(corr_replicates))
df_sameanchor_diffcellline <- data.frame(List = "Same anchor\ndifferent cell lines", Correlation = unlist(corr_sameanchor_diffcellline))
df_diffanchors_samecellline <- data.frame(List = "Different anchors\nsame cell line", Correlation = unlist(corr_diffanchors_samecellline))
df_diffanchors_diffcellline <- data.frame(List = "Different anchors\ndifferent cell lines", Correlation = unlist(corr_diffanchors_diffcellline))
df_meananchors <- data.frame(List = "Different anchors\nacross cell lines", Correlation = unlist(corr_meananchors))
df_meancelllines <- data.frame(List = "Across anchors\ndifferent cell lines", Correlation = unlist(corr_meancelllines))

mean(df_meancelllines$Correlation)

corr_data <- rbind(df_replicates, df_diffanchors_diffcellline, df_sameanchor_diffcellline, df_diffanchors_samecellline, df_meananchors, df_meancelllines)
corr_data$List <- factor(corr_data$List, levels = rev(c("Replicates", "Different anchors\ndifferent cell lines", "Same anchor\ndifferent cell lines", "Different anchors\nsame cell line", "Different anchors\nacross cell lines", "Across anchors\ndifferent cell lines")))

colours_correlations <- c( "#4C5454","#006633","#1D71B8","#82368C","#F39200","#BE1622") 


p <- ggplot(corr_data, aes(x = Correlation, y = List, color = List)) +
  geom_jitter(size = 1, width = 0, height = 0.2, show.legend = FALSE) +
  scale_color_manual(values = colours_correlations) +
  geom_boxplot(alpha = 0.3, color = "black", width = 0.6, outlier.shape = NA, linewidth = 0.1) +
  labs(x = "Correlation", y = "List", color = "List", fill = "List") +
  xlim(0, 1) +
  theme_pie

plot_correlation_2D <- function(data, x_var, y_var, x_label, y_label, x_range = c(-10, 10), y_range = c(-10, 10)) {
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_bin2d(bins = 150, show.legend = FALSE) +
    scale_fill_continuous(type = "viridis") +
    coord_cartesian(xlim = x_range, ylim = y_range) +
    labs(x = x_label, y = y_label, fill = "Density") +
    annotate("text", x = x_range[2] - 1.5, y = y_range[1] + 0.5, label = paste("R = ", round(cor(data[[x_var]], data[[y_var]], method = "pearson", use = "complete.obs"), 2)), family = "Helvetica", size = 2) +
    theme_pie}

p <- plot_correlation_2D(data = dipper_all,
                         x_var = "L2FC_HAP", y_var = "L2FC_HEK",
                         x_label = "L2FC - HAP1", y_label = "L2FC - HEK293T",
                         x_range = c(-10, 6), y_range = c(-10, 6))


####### Short deletions
arrayed_short <- read_csv("/Users/juliane/Documents/dels/files/FACS/ShortDeletions_statistics.csv") %>%
  mutate(forward_pegRNA = factor(forward_pegRNA, levels = c("DPf06", "DPf07", "DPf12", "DPf09", NA)))

p <- ggplot(arrayed_short %>% filter(Day == "D7"), aes(x = reverse_pegRNA, y = forward_pegRNA)) +
  geom_tile(aes(fill = DeletionRate_mean), color = "white") +
  scale_fill_gradient(low = "#d3d3d3", high = "#b2182b") +
  geom_text(aes(label = paste0(round(DeletionRate_mean, 0), "%")), vjust = 1, size = 5/.pt) +  # Add annotations
  labs(x = "Reverse pegRNA", y = "Forward pegRNA", fill = "Deletion rate") +
  theme_pie +
  #theme(legend.position = "right") +
  theme(axis.text.y = element_text(size = 7, hjust = 0),  # Adjust size of y-axis labels
        axis.title.x = element_blank(),                   # Remove x-axis title
        axis.title.y = element_text(size = 7, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        aspect.ratio = 1)  # Ensure squares remain square

####### InDel rate

##### Import data
indel <- read_csv("supp_indel_statistics.csv")

# Prepare data for the stacked and individual bars
tidy_data <- indel %>%
  pivot_longer(cols = starts_with("perc"), names_to = "category", values_to = "percentage") %>%
  mutate(sample = case_when(str_detect(category, "DPf25") ~ "DPf25_DPr08", str_detect(category, "DPf27") ~ "DPf27_DPr08"),
    type = if_else(groupname == "correct", "correct", "others"))

# Prepare data for independent samples
independent_data <- tidy_data %>%
  mutate( groupname = case_when(groupname %in% c("dnt_mutated", "abnormal_attB_length") ~ "mutated_attB", # I decided to merge them as there is not much variance
      TRUE ~ groupname),
      groupname = case_when(groupname %in% c("targetsites_mutated", "targetsites_mutated_attBpresent") ~ "targetsites_mutated", # I decided to merge them as there is not much variance
      TRUE ~ groupname),
    sample_group = case_when(
      type == "correct" & sample == "DPf25_DPr08" ~ "DPf25_correct", type == "correct" & sample == "DPf27_DPr08" ~ "DPf27_correct",
      type == "others" & sample == "DPf25_DPr08" ~ "DPf25_others", type == "others" & sample == "DPf27_DPr08" ~ "DPf27_others")) %>%
  group_by(sample_group, groupname) %>%
  summarise(percentage = sum(percentage), .groups = "drop")

# Plotting
p <- ggplot(independent_data, aes(x = sample_group, y = percentage, fill = groupname)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, show.legend = FALSE, color = "black", linewidth = 0.25) +
  #geom_text(aes(label = round(percentage, 1),group = groupname),position = position_stack(vjust = 0.5),size = 2, color = "black") +
  labs(x = "Sample", y = "Percentage", fill = "Outcome") +
  scale_fill_manual(values = c("correct" = "#d3d3d3", "mutated_attB" ="#BE1622","targetsites_mutated" = "#F9B233")) +
  theme_pie

##########################################################################
################################ FIGURE 2 ################################
####################  determinants of editing  ###########################
##########################################################################

## Build the model in python for nice plotting of SHAP (DP9)

####### Characterize main determinants ----

#### Length plot
dipper <- dipper %>%
  mutate(start_mS_kb = start_mS / 1000)

p <- ggplot(dipper %>% filter(length_mS_DPr04 > 0, length_mS_DPr04 < 10000), aes(x = -start_mS_kb, y = L2FC_DPr04)) +
  geom_point(shape = 16, alpha = 0.8, size = 0.2, show.legend = FALSE, color = "#1D71B8") +
  labs(x = "Deletion position on Chr13 [kb] (reverse)", y = "Log2 fold change for DPr04") +
  scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5),labels = c( "-5", "-2.5", "0", "2.5", "5", "7.5"), limits = c(-5, 7.5), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(-99825, -99822.5, -99820, -99817.5, -99815),labels = c( "99825", "", "99820", "", "99815"), limits = c(-99826, -99814.5), expand = c(0, 0)) +
  theme_pie

# also plot the benchling annotations
regions <- data.frame(Region = c("mScarlet", "CAGG", "CLYBL homology"), Start = c(-99824.702, -99822.961, -99822.145), End = c(-99825.477, -99824.702, -99822.961))
p <- ggplot(regions, aes(y = Region)) +
  geom_segment(aes(x = Start, xend = End, yend = Region), size = 3, color = "blue") +
  scale_x_continuous(breaks = c(-99825, -99822.5, -99820, -99817.5, -99815),labels = c( "99825", "", "99820", "", "99815"), limits = c(-99826, -99814.5), expand = c(0, 0)) +
  theme_pie

#### Length plot by category
# mead editing rates per length bin
data_bins <- dipper %>%
  mutate(length_bin2 = cut(length_mS_DPr04, breaks = c(-0.1, 0 ,1000, 4999, 9999, 99999, 499999, 999999, Inf), labels = c("neg ctrl", "<1", "1-5", "5-10", "10-100", "100-500", "500-1000", ">1000")),
         length_bin3 = cut(length_mS_DPr08, breaks = c(-0.1, 0 ,1000, 4999, 9999, 99999, 499999, 999999, Inf), labels = c("neg ctrl", "<1", "1-5", "5-10", "10-100", "100-500", "500-1000", ">1000")),
         length_bin4 = cut(length_mS_BGHr01, breaks = c(-0.1, 0 ,1000, 4999, 9999, 99999, 499999, 999999, Inf), labels = c("neg ctrl", "<1", "1-5", "5-10", "10-100", "100-500", "500-1000", ">1000")))

data_bins %>%
  group_by(length_bin2) %>%
  summarise(mean = mean(L2FC,na.rm = TRUE))

data_DPr04 <- data_bins %>%
  dplyr::select(length_bin = length_bin2, L2FC = L2FC_DPr04) %>%
  mutate(Condition = "DPr04")
data_DPr08 <- data_bins %>%
  dplyr::select(length_bin = length_bin3, L2FC = L2FC_DPr08) %>%
  mutate(Condition = "DPr08")
data_BGHr01 <- data_bins %>%
  dplyr::select(length_bin = length_bin4, L2FC = L2FC_BGHr01) %>%
  mutate(Condition = "BGHr01")
data_long <- bind_rows(data_DPr04, data_DPr08, data_BGHr01) %>%
  mutate(Condition = factor(Condition, levels = c("DPr04", "DPr08", "BGHr01")))

data_long_selected <- data_long %>% filter(Condition == "DPr04")
anova_result <- aov(L2FC ~ length_bin, data = data_long_selected)
summary(anova_result)
summary(aov(L2FC ~ length_bin, data = data_long_selected))

tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
plot(tukey_result)

custom_colors <- c("DPr04"  = "#A0D2F3", "DPr08"  = "#AC8BBF", "BGHr01" = "#7BA68F")
neg_control_medians <- data_long %>%
  filter(length_bin == "neg ctrl") %>%
  group_by(Condition) %>%
  summarise(median_L2FC = median(L2FC, na.rm = TRUE))

p <- ggplot(data_long, aes(x = as.factor(length_bin), y = L2FC, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.1, position = position_dodge(width = 0.8), show.legend = FALSE) +
  labs(x = "Length Category [kb]", y = "L2FC", fill = "Condition") +
  geom_hline(data = neg_control_medians, aes(yintercept = median_L2FC, color = Condition), linetype = "dashed", size = 0.25, show.legend = FALSE) +
  scale_fill_manual(values = custom_colors) +
  coord_cartesian(xlim = c(0.5, 8.5), ylim = c(-5, 7.5), expand = FALSE) +
  theme_pie


#### Correlation with prediction tools
dset <- dipper %>%
  filter(length_mS_DPr04 < 10000) %>%
  filter(targetgroup != "control")

p <- ggplot(dset, aes(x = L2FC, y = pridict2_score_mS_HEK)) +
  geom_smooth(method = "lm", se = FALSE, color = "#D5CCAE", size = 0.2) +
  #geom_bin2d(bins = 100, show.legend = FALSE) +
  geom_point(shape = 16, alpha = 0.6, size = 0.8) +
  scale_fill_continuous(type = "viridis") +
  coord_cartesian(xlim = c(-4,7.5), ylim =  c(0, 40), expand = FALSE) +
  labs(x = "L2FC for deletions < 10kb", y = "PRIDICT2 score", fill = "Density") +
  annotate("text", x = 6, y = 1, label = paste("R = ", round(cor(dset$L2FC, dset$pridict2_score_mS_HEK, method = "pearson", use = "complete.obs"), 2)), family = "Helvetica", size = 2) +
  theme_pie
pdf(file = paste0(outdir, "Dipper_results_PRIDICT2.pdf"), width = 3 / 2.54, height = 3 / 2.54); grid.draw(set_panel_size(p, width = unit(3, "cm"), height = unit(3, "cm"))); dev.off()

p <- ggplot(dset, aes(x = L2FC, y = percGC_pbs)) +
  geom_jitter(width = 0.01, shape = 16, size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "#D5CCAE", size = 0.5) +
  labs(x = "L2FC", y = "%GC (PBS)") +
  coord_cartesian(ylim = c(0,1), xlim =  c(-2.5, 7.5), expand = FALSE) +
  annotate("text", x = 5, y = 0.125, label = paste("R = ", round(cor(dset$L2FC, dset$percGC_pbs, method = "pearson", use = "complete.obs"), 2)), family = "Helvetica", size = 2) +
  theme_pie
pdf(file = paste0(outdir, "Dipper_results_percGC_PBS.pdf"), width = 3 / 2.54, height = 3 / 2.54); grid.draw(set_panel_size(p, width = unit(3, "cm"), height = unit(3, "cm"))); dev.off()


p <- ggplot(dset, aes(x = L2FC, y = DeepCas9)) +
  geom_smooth(method = "lm", se = FALSE, color = "#D5CCAE", size = 0.2) +
  #geom_bin2d(bins = 100, show.legend = FALSE) +
  geom_point(shape = 16, alpha = 0.6, size = 0.8) +
  scale_fill_continuous(type = "viridis") +
  coord_cartesian(xlim = c(-4,7.5), ylim =  c(0, 80), expand = FALSE) +
  labs(x = "L2FC for deletions < 10kb", y = "DeepCas9 score", fill = "Density") +
  annotate("text", x = 6, y = 2, label = paste("R = ", round(cor(dset$L2FC, dset$DeepCas9, method = "pearson", use = "complete.obs"), 2)), family = "Helvetica", size = 2) +
  theme_pie
pdf(file = paste0(outdir, "Dipper_results_DeepCas9.pdf"), width = 3 / 2.54, height = 3 / 2.54); grid.draw(set_panel_size(p, width = unit(3, "cm"), height = unit(3, "cm"))); dev.off()


#### Contact frequency
dset <- dipper %>%
  filter(length_mS_DPr04 > 100000) %>%
  filter(targetgroup != "control")

p <- ggplot(dset, aes(x = L2FC, y = MCC_mS_coverage_bin1000)) +
  geom_smooth(method = "lm", se = FALSE, color = "#D5CCAE", size = 0.2) +
  #geom_bin2d(bins = 150, show.legend = FALSE) +
  geom_point(shape = 16, alpha = 0.4, size = 0.8) +
  scale_fill_continuous(type = "viridis") +
  coord_cartesian(xlim = c(-5.5,3), ylim =  c(-10, 400), expand = FALSE) +
  labs(x = "L2FC for deletions > 100kb", y = "Contact frequency") +
  annotate("text", x = 2, y = 0, label = paste("R = ", round(cor(dset$L2FC, dset$MCC_mS_coverage_bin1000, method = "pearson", use = "complete.obs"), 2)), family = "Helvetica", size = 2) +
  theme_pie
pdf(file = paste0(outdir, "Dipper_results_MCC.pdf"), width = 3 / 2.54, height = 3 / 2.54); grid.draw(set_panel_size(p, width = unit(3, "cm"), height = unit(3, "cm"))); dev.off()

dset <- dipper %>%
  filter(length_mS_DPr04 > 0) %>%
  filter(targetgroup != "control")

p <- ggplot(dset %>% filter(MCC_mS_coverage_bin1000 > 0), aes(x = L2FC, y = MCC_mS_coverage_bin1000)) +
  geom_smooth(method = "lm", se = FALSE, color = "#D5CCAE", size = 0.2) +
  geom_point(shape = 16, alpha = 0.4, size = 0.8) +
  scale_y_log10() +
  scale_fill_continuous(type = "viridis") +
  #coord_cartesian(xlim = c(-5.5, 3), expand = FALSE) +
  labs(x = "L2FC", y = "Contact frequency (log scale)") +
  annotate("text", x = 2, y = 0, label = paste("R = ", round(cor(dset$L2FC, dset$MCC_mS_coverage_bin1000, method = "pearson", use = "complete.obs"), 2)), family = "Helvetica", size = 2) +
  theme_pie
pdf(file = paste0(outdir, "Dipper_results_MCC_alldata.pdf"), width = 3 / 2.54, height = 3 / 2.54); grid.draw(set_panel_size(p, width = unit(3, "cm"), height = unit(3, "cm"))); dev.off()


#### Plot predicted vs measured
test_data_corefeatures <- read.csv("/Users/juliane/Documents/dels/files/modelling/dipper_test_with_results_core.csv")

p <- ggplot(test_data_corefeatures, aes(x = L2FC, y = L2FC_pred)) +
  #geom_smooth(method = "lm", se = FALSE, color = "#D5CCAE", size = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "#D5CCAE", linetype = "dashed") +
  #geom_bin2d(bins = 75, show.legend = FALSE) +
  geom_point(shape = 16, alpha = 0.4, size = 0.8, color = "black") +
  scale_fill_continuous(type = "viridis") +
  coord_cartesian(xlim = c(-5,7.5), ylim =  c(-5, 7.5), expand = FALSE) +
  labs(x = "Observed L2FC", y = "Predicted L2FC", fill = "Density") +
  annotate("text", x = 6, y = -4.75, label = paste("R = ", round(cor(test_data_corefeatures$L2FC_pred, test_data_corefeatures$L2FC, method = "pearson", use = "complete.obs"), 2)), family = "Helvetica", size = 2) +
  theme_pie
pdf(file = paste0(outdir, "Dipper_results_modeling_testset.pdf"), width = 4 / 2.54, height = 4 / 2.54); grid.draw(set_panel_size(p, width = unit(4, "cm"), height = unit(4, "cm"))); dev.off()

view(test_data_corefeatures %>% select(L2FC, L2FC_pred, length_mS_DPr04, deepprime_score_mS, ATAC_HAP1, MCC_mS_coverage_bin1000, pridict2_score_mSd_HEK) %>% filter(L2FC_pred > 2))


##########################################################################
################################ FIGURE 3 ################################
#################  ultra long deletions and ddPCR  ########################
##########################################################################

### ultra long deletions
guides <- read_csv("/Users/juliane/Library/CloudStorage/GoogleDrive-jw38@sanger.ac.uk/My Drive/PhD/41_PDel/30_Dipper_screen/followup/20230927_LongerDeletions_guides.csv")

# get data for HAP1 DPr08 deletions: this is the prepared file so can skip the next rows
arrayed_repeat <- read_csv("/Users/juliane/Documents/dels/files/FACS/ArrayedLongDeletions_statistics.csv")

##### Prepare data - later only read in result ######

# Import the FACS: long deletions for HAP1 DPr08 (round 1)
arrayed_repeat1 <- read_csv("/Users/juliane/Library/CloudStorage/GoogleDrive-jw38@sanger.ac.uk/My Drive/PhD/41_PDel/rawdata/FACS/20240611_arrayed_deletions_analysis/ANALYSISREADY_Arrayed_short_deletions/Arrayed_deletion_ultralong_DPr08_analysisNov5.csv") %>%
  separate("...1", into = c("temp1","temp2", "WellID", "forward_pegRNA", "reverse_pegRNA"), sep = "_") %>%
  mutate(reverse_pegRNA = gsub(".fcs", "", reverse_pegRNA), Day = temp2, Cellline = ifelse(temp1 == "HAP", "HAP1", "HEK293T")) %>%
  mutate("DeletionRate" = Freq_mS_neg) %>%
  group_by(forward_pegRNA, reverse_pegRNA, Day, Cellline) %>%
  summarise(DeletionRate = mean(DeletionRate, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Repeat = "R1") %>%
  left_join(guides %>% dplyr::select("Name","length_DPr03"), by = c("forward_pegRNA" = "Name")) %>%
  filter(Cellline == "HAP1", Day == "D7")

# Import the FACS: long deletions repeat DPr08 (round 2)
arrayed_repeat2 <- read_csv("/Users/juliane/Library/CloudStorage/GoogleDrive-jw38@sanger.ac.uk/My Drive/PhD/41_PDel/rawdata/FACS/20241031_ArrayedDeletionRepeats/LongDeletions/LongDeletions_statistics.csv") %>%
  separate("...1", into = c("temp1","temp2", "WellID","forward_pegRNA", "reverse_pegRNA", "Replicate"), sep = "_") %>%
  mutate(Replicate = gsub(".fcs", "", Replicate)) %>%
  mutate("DeletionRate" = Freq_mS_neg) %>%
  group_by(forward_pegRNA, reverse_pegRNA, Day, Cellline) %>%
  summarise(DeletionRate = mean(DeletionRate, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Repeat = "R2") %>%
  left_join(guides %>% dplyr::select("Name","length_DPr03"), by = c("forward_pegRNA" = "Name")) %>%
  filter(Cellline == "HAP1", Day == "D7")

# Combine the two DPr08 datasets
arrayed_repeat <- rbind(arrayed_repeat1, arrayed_repeat2)

# Adjust length because it is DPR03
arrayed_repeat$deletion_length <- arrayed_repeat$length_DPr03
arrayed_repeat$deletion_length <- arrayed_repeat$deletion_length +334
arrayed_repeat$deletion_length[arrayed_repeat$forward_pegRNA == "None"] <- 0
arrayed_repeat$deletion_length[arrayed_repeat$reverse_pegRNA == "None"] <- 0
arrayed_repeat$length_label <- round(arrayed_repeat$deletion_length / 1000, 0)

arrayed_repeat <- arrayed_repeat %>% dplyr::select(-length_DPr03)
arrayed_repeat <- arrayed_repeat %>% filter(length_label!= "BFP")

# Create a unique name for each pegRNA pair
arrayed_repeat$pegRNA_pair <- paste(arrayed_repeat$reverse_pegRNA,arrayed_repeat$forward_pegRNA, sep = "_")


####### Plot 

# Calculate the mean deletion rate for each pegRNA pair
arrayed_repeat_grouped <- arrayed_repeat %>%
  group_by(pegRNA_pair) %>%
  summarise(DeletionRate_mean = mean(DeletionRate, na.rm = TRUE), DeletionRate_std = sd(DeletionRate, na.rm = TRUE)) %>%
  ungroup()

# Calculate the significance of controls vs 1Mb+ deletions
controls <- arrayed_repeat %>%
  filter(pegRNA_pair %in% c("None_DPf17", "None_DPf18", "None_None", "DPr08_None")) %>%
  pull(DeletionRate)
treated <- arrayed_repeat %>%
  filter(pegRNA_pair %in% c("DPr08_DPf19", "DPr08_DPf25", "DPr08_DPf40")) %>%
  pull(DeletionRate)
t_test_result <- t.test(treated, controls)
print(t_test_result$p.value)

pegRNA_levels <- arrayed_repeat %>% arrange(deletion_length) %>% pull(pegRNA_pair) %>% unique()
labels_temp <- arrayed_repeat %>% arrange(deletion_length) %>% pull(length_label) %>% unique()

labels <- c("only DPf17", "DPr08", "WT", "only DPf18", labels_temp)
labels <- labels[labels != 0] # fixing them manually

p <- ggplot() +
  geom_bar(data = arrayed_repeat_grouped, 
           aes(x = factor(pegRNA_pair, levels = pegRNA_levels), y = DeletionRate_mean),
           stat = "identity", position = position_dodge(width = 1), width = 0.8) +  # Bar plot with dodge
  geom_point(data = arrayed_repeat, 
             aes(x = factor(pegRNA_pair, levels = pegRNA_levels), y = DeletionRate),
             size = 0.6, alpha = 1, stroke = 0) +
  labs(x = "Deletion length [kb]", y = "Deletion rate [%]") +
  scale_x_discrete(labels = labels) +
  coord_cartesian(ylim = c(0, 80)) +
  theme_pie

### ddPCR
data_Rho <- read_csv("ddPCR_Rho_data.csv")
data_DMD <- read_csv("ddPCR_DMD_data.csv")

#### Plot normalized values
sorting_order = c("Hap1_PE7", "Rho1", "Rho3", "Rho4", "Rho5", "Rho7", "Rho8", "Rho9",  "Rho11")
p <- ggplot(data = data_Rho_mean %>% filter(probe == "Probe Exon1"), aes(x = factor(sample, levels = sorting_order), y = deletionrate_norm)) +
  stat_summary(fun = "mean", geom = "bar", fill = "#A0D2F3", color = "black", width = 0.7) +
  geom_jitter(width = 0.15, size = 1, shape = 21, fill = "white", color = "black") +
  scale_y_continuous(limits = c(-5, 100), breaks = seq(0, 100, by = 10), labels = ~ ifelse(.x %% 20 == 0, as.character(.x), ""), expand = c(0, 0)) +
  labs(x = "Deletion", y = "Normalized deletion rate [%]") +
  theme_pie

p <- ggplot(data = data_Rho_mean %>% filter(probe == "Probe Exon4"), aes(x = factor(sample, levels = sorting_order), y = deletionrate_norm)) +
  stat_summary(fun = "mean", geom = "bar", fill = "#A0D2F3", color = "black", width = 0.7) +
  geom_jitter(width = 0.15, size = 1, shape = 21, fill = "white", color = "black") +
  scale_y_continuous(limits = c(-5, 100), breaks = seq(0, 100, by = 10), labels = ~ ifelse(.x %% 20 == 0, as.character(.x), ""), expand = c(0, 0)) +
  labs(x = "Deletion", y = "Normalized deletion rate [%]") +
  theme_pie

sorting_order = c("Hap1_PE7", "DMD4", "DMD5", "DMD6")
p <- ggplot(data = data_DMD_mean, aes(x = factor(sample, levels = sorting_order), y = deletionrate_norm, fill = probe)) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.7), color = "black", width = 0.7, show.legend = FALSE) +
  geom_jitter(aes(color = probe), position = position_dodge(width = 0.7), size = 1, shape = 21, fill = "white", show.legend = FALSE) +
  scale_y_continuous(limits = c(-5, 100), breaks = seq(0, 100, by = 10), labels = ~ ifelse(.x %% 20 == 0, as.character(.x), ""), expand = c(0, 0)) +
  labs(x = "Deletion", y = "Normalized deletion rate [%]") +
  scale_fill_manual(values = c("#A0D2F3", "#1D71B8")) +
  scale_color_manual(values = c("black", "black")) +
  theme_pie


