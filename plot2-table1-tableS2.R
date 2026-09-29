#===== workspace
# Figure 2 (main paper), Table 1 (main paper) and Table S2 (supporting info):
# cross-region validation, settings E1-E4 (each latitude zone used in turn as
# the test region, the remaining zones used for training).
#   in:  output/RData/analysis_cross_region.RData
#   out: output/figure/plot_lat-abs_diff1.pdf   (Figure 2)
#        Table 1  -> printed to the console as LaTeX
#        Table S2 -> printed to the console as LaTeX
rm(list = ls())

library(ggplot2)
library(tidyverse)
library(xtable)
library(reshape2)

#===== set-up
load("output/RData/analysis_cross_region.RData")

# E1 = Tropics, E2 = Subtropics, E3 = Mid-Latitudes, E4 = Polar (as test region)
zone_labels <- c("Tropics", "Subtropics", "Mid-Lat", "Polar")
Settings <- paste0("E", 1:nrow(eval))

#===========================================================
# Table 1: sample sizes and split configuration for E1-E4
#===========================================================
data_plot_aux <- as.data.frame(cbind(Settings, round(eval[, c("n_train", "n_test", "lat_l", "lat_u")], 2)))
data_plot_aux_m <- melt(data_plot_aux, id.vars = "Settings")     # long format
data_plot_aux_w <- dcast(data_plot_aux_m, variable ~ Settings)   # wide format (transpose style)

cat("\n===== Table 1 (main paper): sample sizes, cross-region validation =====\n")
print(xtable(data_plot_aux_w, caption = "", label = ""), include.rownames = FALSE)

#===========================================================
# Table S2: delta MSE (SE, SE_sub) for LM and 1NN, E1-E4
#===========================================================
method_plot       <- c("lm_mse", "nn_mse")
method_plot_short <- c("LM", "1NN")
method_se         <- c("lm_se", "nn_se")
method_sub_se     <- c("lm_sub_se", "nn_sub_se")

table_s2 <- data.frame(Settings)

## 1. Method columns: estimate (SE_full, SE_sub)
for (i in seq_along(method_plot)) {
  table_s2[[method_plot_short[i]]] <-
    sprintf("%.3f (%.3f, %.3f)",
            as.numeric(eval_diff[, method_plot[i]]),
            as.numeric(eval_diff[, method_se[i]]),      # SE on full test
            as.numeric(eval_diff[, method_sub_se[i]]))  # SE on subsample
}

## 2. Add sample size columns once per row
table_s2$n_test    <- eval_diff[, "n_test"]
table_s2$n_subtest <- eval_diff[, "n_sub_test"]

## 3. Reorder columns
table_s2 <- table_s2[, c("Settings", "n_test", "n_subtest", method_plot_short)]

cat("\n===== Table S2 (supporting info): delta MSE, cross-region validation =====\n")
print(xtable(table_s2, caption = "", label = ""), include.rownames = FALSE)

#===========================================================
# Figure 2: MSE differences vs OC4 with 95% CI, one facet per test zone
#===========================================================
eval_diff <- as.data.frame(eval_diff)
eval_diff$Settings <- paste0(zone_labels, " test")

target_settings <- paste0(zone_labels, " test")

use_se_type <- "sub_se"   # conservative subsample-based SE (see Section S3)
zval <- 1.96
method_key <- c(lm = "LM", nn = "1NN")

df_long <- as_tibble(eval_diff) %>%
  mutate(Settings = factor(Settings, levels = target_settings)) %>%
  pivot_longer(
    cols = matches("^(lm|nn)_(mse|se|sub_se)$"),
    names_to = c("m", ".value"),
    names_pattern = "^(lm|nn)_(mse|se|sub_se)$"
  ) %>%
  mutate(
    Method = recode(m, !!!method_key),
    est    = as.numeric(mse),
    se_use = as.numeric(.data[[use_se_type]]),
    lo     = est - zval * se_use,
    hi     = est + zval * se_use,
    Method = factor(Method, levels = c("LM", "1NN"))
  )

col_point <- "grey55"
col_bar   <- "grey55"
col_dash  <- "grey30"

p <- ggplot(df_long, aes(x = Method, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9, colour = col_dash) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = 0.28, linewidth = 1.1, colour = col_bar) +
  geom_point(shape = 16, size = 3.2, colour = col_point) +
  facet_wrap(~ Settings, ncol = 4) +
  xlab("") +
  ylab("MSE Difference Compared to OC4 \n MSE(Method) - MSE(OC4)") +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 15, colour = "black", face = "bold"),
    axis.title.y = element_text(size = 13),
    axis.text.x  = element_text(size = 11, colour = "black"),
    axis.text.y  = element_text(size = 11)
  )

save_file <- "output/figure/plot_lat-abs_diff1.pdf"
ggsave(save_file, plot = p, width = 8, height = 4)
