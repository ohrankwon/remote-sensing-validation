#===== workspace
# Figure 1 (main paper): within-region validation, mean percent change in MSE
# for each algorithm when using cluster-based splits (Clust) instead of
# observation-level random splits (NoClust).
#   in:  output/RData/analysis_random_mse_clust.RData              (Clust)
#        output/RData/analysis_random_mse_noclust.RData (NoClust)
#   out: output/figure/plot_random_80_pctchange_mean_se.pdf
rm(list = ls())

library(ggplot2)
library(tidyverse)
library(scales)

#===============================
# 1. Common setup
#===============================
which_split <- "R"

method_plot <- c("lm_mse", "oc4_mse", "nn_mse")
new_method_names <- c("LM", "OC4", "1NN")
name_map <- setNames(new_method_names, method_plot)

#===============================
# 2. WITH ST-cluster (Clust)
#===============================
load("output/RData/analysis_random_mse_clust.RData")

Settings <- paste0(which_split, 1:nrow(eval))
data_plot <- as.data.frame(cbind(Settings = Settings, eval[, method_plot]))

#===============================
# 3. WITHOUT ST-cluster (NoClust)
#===============================
load("output/RData/analysis_random_mse_noclust.RData")

Settings2 <- paste0(which_split, 1:nrow(eval))
data_plot2 <- as.data.frame(cbind(Settings = Settings2, eval[, method_plot]))

#===============================
# 4. Percent change (Clust vs NoClust)
#===============================
num_cols <- method_plot

data_plot[num_cols]  <- lapply(data_plot[num_cols],  as.numeric)
data_plot2[num_cols] <- lapply(data_plot2[num_cols], as.numeric)

# Compute percent change: (with - without)/without * 100
data_plot3 <- data_plot2
data_plot3[num_cols] <- (data_plot[num_cols] - data_plot2[num_cols]) /
  data_plot2[num_cols] * 100

data_plot3_long <- data_plot3 %>%
  pivot_longer(
    cols = -Settings,
    names_to = "Method",
    values_to = "PctChange"
  )

data_plot3_long$Method    <- name_map[data_plot3_long$Method]
data_plot3_long$PctChange <- as.numeric(data_plot3_long$PctChange)

#===============================
# 5. Percent change: mean and 95% CI by method
#     The interval is Monte Carlo error across the num_repeat split repeats.
#     It is NOT the subsample SE used in Figures 2 and S5, which measures
#     sampling error across test observations within a single split.
#===============================
data_plot3_long <- data_plot3_long %>%
  mutate(
    Method = factor(Method, levels = c("OC4", "LM", "1NN"))
  )

summary_pct <- data_plot3_long %>%
  group_by(Method) %>%
  summarise(
    mean_pct = mean(PctChange, na.rm = TRUE),
    sd_pct   = sd(PctChange,   na.rm = TRUE),
    n        = sum(!is.na(PctChange)),
    se_pct   = sd_pct / sqrt(n),
    ci_pct   = 1.96 * se_pct,
    .groups  = "drop"
  )

#===============================
# 6. Plot: mean percent change with 95% CI
#===============================
p_pct <- ggplot(summary_pct,
                aes(x = Method, y = mean_pct)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) +
  geom_point(size = 3, stroke = 1.2, colour = "grey50") +
  geom_errorbar(aes(ymin = mean_pct - ci_pct,
                    ymax = mean_pct + ci_pct),
                width = 0.25,
                linewidth = 1.1,
                colour = "grey50") +
  ylab("Mean percent change in MSE\n(Clust - NoClust) / NoClust * 100") +
  xlab("") +
  scale_y_continuous(labels = function(x) sprintf("%.1f%%", x)) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title.y    = element_text(size = 13),
    axis.text.x     = element_text(size = 12),
    axis.text.y     = element_text(size = 12)
  )

ggsave("output/figure/plot_random_80_pctchange_mean_se.pdf",
       plot = p_pct, width = 4, height = 4)

print(summary_pct)
