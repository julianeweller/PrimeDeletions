# Scripts for thesis plots
library(tidyverse)
library(dplyr)
library(ggplot2)
library(egg)
library(grid)
library(scales)
library(reshape2)

path = "/"
outdir = paste0(path, "/figures/")
outdir_supp = paste0(path, "/supplements/")

theme_pie <-   theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        text = element_text(family = "Helvetica"),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = 8),
        plot.background = element_blank())

################# Load data #################

data <- read_csv("data.csv") 

################# Figure 4 #################

####### 4d) replicate correlation -----
plot_correlation_2D <- function(data, x_var, y_var, x_label, y_label, x_range = c(-10, 5), y_range = c(-10, 5)) {
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_bin2d(bins = 250, show.legend = FALSE) +
    scale_fill_continuous(type = "viridis") +
    coord_cartesian(xlim = x_range, ylim = y_range) +
    labs(x = x_label, y = y_label, fill = "Density") +
    annotate("text", x = x_range[2] - 1.5, y = y_range[1] + 0.5, label = paste("R = ", round(cor(data[[x_var]], data[[y_var]], method = "pearson", use = "complete.obs"), 2)), family = "Helvetica", size = 2) +
    theme_pie}

p <- plot_correlation_2D(data = data,
                         x_var = "L2FC_R1", y_var = "L2FC_R2",
                         x_label = "Replicate 1 [L2FC]", y_label = "Replicate 2 [L2FC]")


####### 4e) replicate correlation between different subsets -----
corr_replicates <- list(
  replicates_D11 = cor(data$`L2FC_1-D11D4`, data$`L2FC_2-D11D4`, use = "complete.obs", method = "pearson"),
  replicates_D14 = cor(data$`L2FC_1-D14D4`, data$`L2FC_2-D14D4`, use = "complete.obs", method = "pearson"),
  replicates_D21 = cor(data$`L2FC_1-D21D4`, data$`L2FC_2-D21D4`, use = "complete.obs", method = "pearson"),
  replicates_D28 = cor(data$`L2FC_1-D28D4`, data$`L2FC_2-D28D4`, use = "complete.obs", method = "pearson"))
corr_days <- list(
  days_D11_D14 = cor(data$L2FC_D11D4, data$L2FC_D14D4, use = "complete.obs", method = "pearson"),
  days_D11_D21 = cor(data$L2FC_D11D4, data$L2FC_D21D4, use = "complete.obs", method = "pearson"),
  days_D11_D28 = cor(data$L2FC_D11D4, data$L2FC_D28D4, use = "complete.obs", method = "pearson"),
  days_D14_D21 = cor(data$L2FC_D14D4, data$L2FC_D21D4, use = "complete.obs", method = "pearson"),
  days_D14_D28 = cor(data$L2FC_D14D4, data$L2FC_D28D4, use = "complete.obs", method = "pearson"),
  days_D21_D28 = cor(data$L2FC_D21D4, data$L2FC_D28D4, use = "complete.obs", method = "pearson"))
corr_merged_acrossdays <- list(
  merged_acrossdays = cor(data$L2FC_R1, data$L2FC_R2, use = "complete.obs", method = "pearson"))

df_replicates <- data.frame(List = "Replicates\nat each timepoint", Correlation = unlist(corr_replicates))
df_days <- data.frame(List = "Between\ntimepoints", Correlation = unlist(corr_days))
df_merged <- data.frame(List = "Replicates\nacross timepoints", Correlation = unlist(corr_merged_acrossdays))
corr_data <- rbind(df_replicates, df_days, df_merged)
corr_data$List <- factor(corr_data$List, levels = rev(c("Between\ntimepoints", "Replicates\nat each timepoint", "Replicates\nacross timepoints")))

write_csv(corr_data, "dunnock_replicatecorrelations.csv")

colours_correlations <- c( "#4C5454","#006633","#1D71B8") 

p <- ggplot(corr_data, aes(x = Correlation, y = List, color = List)) +
  geom_jitter(size = 1, width = 0, height = 0.15, show.legend = FALSE) +
  scale_color_manual(values = colours_correlations) +
  geom_boxplot(alpha = 0.3, color = "black", width = 0.6, outlier.shape = NA, linewidth = 0.1) +
  labs(x = "Correlation", y = "List", color = "List", fill = "List") +
  xlim(0, 1) +
  theme_pie


####### 4f) Control distributions -----

data_temp <- data %>%
  mutate(group = ifelse(gene_group_label == "NTctrl", "NTctrl", ifelse(gene_group_label == "control", "control", ifelse((positive_control_gene == TRUE &  HAP1_exonic == TRUE), "poscontrol", "main")))) %>%
  mutate(group = factor(group, levels = c("poscontrol","main","control","NTctrl"))) %>%
  select(L2FC, group, dropout, enriched, effected)

p <- ggplot(data_temp, aes(x = L2FC, color = group, fill = group)) + 
  geom_vline(xintercept = 1, linetype = "dotted", color = "black", linewidth = 0.25) + 
  geom_vline(xintercept = -1, linetype = "dotted", color = "black", linewidth = 0.25) + 
  geom_density(aes(y = ..scaled..), alpha = 0.1, trim = FALSE, linewidth = 0.25) + 
  labs(x = "L2FC", y = "Density") + 
  scale_color_manual(values = c("NTctrl" = "darkgray", "control" = "#1D71B8", "poscontrol" = "#006633", "main" = "#BE1622")) +
  scale_fill_manual(values = c("NTctrl" = "white", "control" = "white", "poscontrol" = "#006633", "main" = "#BE1622")) +
  theme_pie + 
  theme(legend.position = c(0.1, 0.9), family = "Helvetica", size = 2)

# count number of oligos per group
table(data_temp$group)

# percentage neither depleted nor enriched
table(data_temp$group, data_temp$effected)[,1] / rowSums(table(data_temp$group, data_temp$effected)) * 100
# percentage depleted
table(data_temp$group, data_temp$dropout)[,2] / rowSums(table(data_temp$group, data_temp$dropout)) * 100
# percentage enriched
table(data_temp$group, data_temp$enriched)[,2] / rowSums(table(data_temp$group, data_temp$enriched)) * 100
# percentage depleted or enriched


####### 4g) Enrichment and depletion analysis -----

p <- ggplot(data, aes(x = rank1, y = L2FC)) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "grey", size = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey", size = 0.2) +
  geom_point(aes(color = colourname), size = 0.5, stroke = 0) +
  scale_color_identity() +
  labs(x = "Rank of pairs\nin dataset", y = "L2FC") +
  coord_cartesian(ylim = c(-7.75, 2.75)) +
  theme_pie
pdf(file = paste0(outdir, "Dunnock_ranking_all.pdf"), width = 3.5 / 2.54, height = 3.5 / 2.54); grid.draw(set_panel_size(p, width = unit(3.5, "cm"), height = unit(3.5, "cm"))); dev.off()

# count and percentage of depleted and enriched
n_depleted <- sum(data$dropout, na.rm = TRUE) # n = 1196
n_enriched <- sum(data$enriched, na.rm = TRUE) # n = 82
n_depleted / nrow(data) * 100
n_enriched / nrow(data) * 100


################# Figure 5 #################

###### 5a: Essentiatility vs dropout

dropout_by_gene <- data %>%
  group_by(main_gene) %>%
  #filter(gencode_coding_annotation == "coding_sequence") %>%
  filter(HAP1_exonic == TRUE) %>%
  summarise(
    dropout_rate = sum(dropout, na.rm = TRUE) / sum(!is.na(dropout)),
    enriched_rate = sum(enriched, na.rm = TRUE) / sum(!is.na(enriched)),
    effected_rate = sum(effected, na.rm = TRUE) / sum(!is.na(effected)),
    n = n()) %>%
  filter( n > 2) %>%
  ungroup() %>%
  arrange(desc(dropout_rate)) %>%
  left_join(data %>% select(main_gene, main_gene_essentiality) %>% distinct(), by = "main_gene")

p <- ggplot(dropout_by_gene, aes(x = main_gene_essentiality, y = dropout_rate)) +
  geom_point(aes(size = n), color = "#BE1622", alpha = 0.3, show.legend = FALSE, stroke = 0) +
  scale_size_continuous(range = c(0.01, 2)) + 
  labs(x = "Essentiality [CRISPR Cas9 knockout]", y = "Dropout rate by gene for deletions targeting coding sequences") +
  theme_pie 

############ 5b: Dropout for control genes

p <- ggplot(data %>% filter(positive_control_gene == TRUE,  HAP1_exonic == TRUE), aes(x = L2FC, y = main_gene)) +
  geom_vline(xintercept = -1, linetype = "dotted", color = "black", linewidth = 0.25) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black", linewidth = 0.25) +
  geom_point(aes(color = colourname), size = 0.5, stroke = 0) +
  scale_color_identity() +
  labs(x = "L2FC", y = "") +
  theme_pie +
  theme(axis.text.y = element_text(size = 6, family = "Helvetica"))

########### 5c: Hit rate for control genes

# plot average hit rate per gene
hit_rate <- data %>% 
  filter(positive_control_gene == TRUE,  HAP1_exonic == TRUE) %>%
  group_by(main_gene) %>%
  summarise(
    dropout_rate = sum(L2FC < -1, na.rm = TRUE) / sum(!is.na(L2FC)),
    n = n()) %>%
  ungroup() %>%
  filter(!is.na(main_gene))

# make bar chart
mean_dropout <- mean(hit_rate$dropout_rate)
median(hit_rate$dropout_rate)
mean_dropout


p <- ggplot(hit_rate, aes(x = main_gene, y = dropout_rate)) +
  geom_hline(yintercept = mean_dropout, linetype = "dotted", color = "black", size = 0.25) +
  geom_bar(stat = "identity", fill = "#F4E8D5") +
  labs(x = "Essential genes", y = "Dropout rate") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_pie +
  coord_flip()

