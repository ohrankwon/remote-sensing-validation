#===== workspace
# Figure S2, Table S1 and Table S3 (supporting information): temporal
# validation, settings T1-T7 (clusters ordered by mean collection date, the
# oldest 20%-80% used for training and the remainder for testing).
#   in:  output/RData/analysis_temporal.RData
#   out: output/figure/plot_time_7.pdf   (Figure S2)
#        Table S1 -> printed to the console as LaTeX
#        Table S3 -> printed to the console as LaTeX
rm(list = ls())

library(ggplot2)
library(tidyverse)
library(scales)
library(xtable)
library(ggpattern)
library(reshape2)

#===== set-up
load("output/RData/analysis_temporal.RData")

which_split <- "T"
Settings <- paste0(which_split, 1:nrow(eval))

# `eval` is a character matrix because lat_l / lat_u hold the split dates
eval2 <- apply(eval, 2, as.numeric)

method_plot       <- c("lm_mse", "oc4_mse", "nn_mse")
method_plot_short <- c("LM", "OC4", "1NN")
method_se         <- c("lm_se", "oc4_se", "nn_se")
new_method_names  <- c("LM", "OC4", "1NN")

method_colors <- c(
  "LM"  = "#7570b3",  # purple
  "OC4" = "#d95f02",  # orange
  "1NN" = "#1b9e77"   # green-teal
)

method_patterns <- c(
  "LM"  = "crosshatch",
  "OC4" = "circle",
  "1NN" = "none"
)

#===========================================================
# Figure S2: MSE by setting and method
#===========================================================
data_plot <- as.data.frame( cbind(Settings, eval[, method_plot] ) )

data_plot_long <- data_plot %>%
  pivot_longer(cols = -Settings,
               names_to = "Method",
               values_to = "Performance")

data_plot_long$Performance <- as.numeric(data_plot_long$Performance)

name_map <- setNames(new_method_names, method_plot)
data_plot_long$Method <- name_map[data_plot_long$Method]

# bars above the y-axis cap are clipped and annotated with their true value
data_plot_long$label <- NA
for(i in 1:nrow(data_plot_long)) {
  if (data_plot_long$Performance[i] >= 0.17) {
    data_plot_long$label[i] <- round(data_plot_long$Performance[i], 2)
    data_plot_long$Performance[i] <- 0.17
  }
}

data_plot_long$Settings <- factor(
  data_plot_long$Settings,
  levels = paste0(which_split, 1:10)
)

data_plot_long$Method <- factor(
  data_plot_long$Method,
  levels = c("LM", "OC4", "1NN")
)

y_min <- 0
y_max <- max(data_plot_long$Performance)

p <- ggplot(
  data_plot_long,
  aes(x = Settings, y = Performance,
      fill = Method, pattern = Method)
) +
  geom_bar_pattern(
    stat     = "identity",
    position = position_dodge(width = 0.8),
    width    = 0.6,
    pattern_fill    = NA,
    pattern_colour  = "grey20",
    pattern_density = 0.3,
    pattern_spacing = 0.06,
    pattern_size    = 0.25,
    colour          = "black",
    size            = 0.2
  ) +
  scale_fill_manual(values  = method_colors,  name = "Method") +
  scale_pattern_manual(values = method_patterns, name = "Method") +
  labs(x = "Settings", y = "MSE") +
  coord_cartesian(ylim = c(y_min, y_max)) +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 0, hjust = 1),
    axis.title      = element_text(size = 14),
    axis.text       = element_text(size = 12),
    legend.title    = element_text(size = 13),
    legend.text     = element_text(size = 12),
    plot.title      = element_text(size = 14, hjust = 0.5)
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  geom_text(
    aes(label = label),
    vjust = 1.2, hjust = -0.3,
    color = "red", fontface = "bold"
  )

save_file <- "output/figure/plot_time_7.pdf"
ggsave(save_file, plot = p, width = 6, height = 4)

#===========================================================
# Table S1: sample sizes and split dates for T1-T7
#   lat_l = last training date, lat_u = first test date
#===========================================================
data_plot_aux <- as.data.frame(cbind(Settings, eval[, c("n_train", "n_test", "lat_l", "lat_u")]))
data_plot_aux_m <- melt(data_plot_aux, id.vars = "Settings")     # long format
data_plot_aux_w <- dcast(data_plot_aux_m, variable ~ Settings)   # wide format (transpose style)

cat("\n===== Table S1 (supporting info): sample sizes, temporal validation =====\n")
print(xtable(data_plot_aux_w, caption = "", label = ""), include.rownames = FALSE)

#===========================================================
# Table S3: delta MSE (SE, SE_sub) for LM and 1NN, T1-T7
#===========================================================
method_plot       <- c("lm_mse", "nn_mse")
method_plot_short <- c("LM", "1NN")
method_se         <- c("lm_se", "nn_se")
method_sub_se     <- c("lm_sub_se", "nn_sub_se")

table_s3 <- data.frame(Settings)

## 1. Method columns: estimate (SE_full, SE_sub)
for (i in seq_along(method_plot)) {
  table_s3[[method_plot_short[i]]] <-
    sprintf("%.3f (%.3f, %.3f)",
            as.numeric(eval_diff[, method_plot[i]]),
            as.numeric(eval_diff[, method_se[i]]),      # SE on full test
            as.numeric(eval_diff[, method_sub_se[i]]))  # SE on subsample
}

## 2. Add sample size columns once per row
table_s3$n_test    <- eval_diff[, "n_test"]
table_s3$n_subtest <- eval_diff[, "n_sub_test"]

## 3. Reorder columns
table_s3 <- table_s3[, c("Settings", "n_test", "n_subtest", method_plot_short)]

cat("\n===== Table S3 (supporting info): delta MSE, temporal validation =====\n")
print(xtable(table_s3, caption = "", label = ""), include.rownames = FALSE)
