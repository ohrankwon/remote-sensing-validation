#===== workspace
# Figure S5 (supporting information): MSE comparison under the combined
# cross-region and temporal generalization scenarios M1 and M2.
#   in:  output/RData/analysis_combined.RData
#   out: output/figure/plot_mode_diff.pdf
rm(list = ls())

library(ggplot2)
library(tidyverse)

#===== set-up
load("output/RData/analysis_combined.RData")

Settings <- paste0("M", 1:nrow(eval))

#========================
# Manual switches
#========================
use_se_type <- "sub_se"   # choose: "se" (full test)  OR  "sub_se" (subsample)
zval <- 1.96

method_key <- c(lm = "LM", nn = "1NN")

#========================
# Build long data
#========================
df_long <- as_tibble(eval_diff) %>%
  mutate(Settings = as.character(Settings)) %>%
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

#========================
# Faceted plot: one facet per setting
#========================
p <- ggplot(df_long, aes(x = Method, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9, colour = col_dash) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = 0.25, linewidth = 1.1, colour = col_bar) +
  geom_point(stroke = 1.2, size = 3, colour = col_point) +
  facet_wrap(~ Settings, ncol = 4) +
  xlab("") +
  ylab("MSE Difference Compared \n to OC4 with 95% CI") +
  theme_bw(base_size = 15) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 15, colour = "black", face = "bold"),
    legend.position = "none",
    axis.title.y = element_text(size = 13),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

save_file <- "output/figure/plot_mode_diff.pdf"
ggsave(save_file, plot = p, width = 4, height = 4)
