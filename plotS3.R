#===== workspace
# Figure S3 (supporting information): cross-region segmented MSE heat map for
# LM and OC4. Each cell is the test MSE for one combination of training (rows)
# and testing (columns) latitude zones; the diagonal is left blank.
#   in:  output/RData/analysis_cross_region_grid.RData
#   out: output/figure/plot_lat-abs-segment_heat-map_16_LM_OC4.pdf
rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)

#===== set-up
load("output/RData/analysis_cross_region_grid.RData")

eval <- as.data.frame(eval)

#-----------------------------
# 1. Attach train/test zone indices to eval
#-----------------------------
pairs <- expand.grid(train_zone = 1:4, test_zone = 1:4)
eval$train_zone <- pairs$train_zone
eval$test_zone  <- pairs$test_zone

zone_labels <- c("Tropics", "Subtropics", "Mid-Lat", "Polar")

#-----------------------------
# 2. Long format, LM and OC4 only
#-----------------------------
heat_long <- eval %>%
  dplyr::select(
    train_zone, test_zone,
    LM  = lm_mse,
    OC4 = oc4_mse
  ) %>%
  tidyr::pivot_longer(
    cols = -(train_zone:test_zone),
    names_to = "Method",
    values_to = "MSE"
  ) %>%
  dplyr::mutate(
    train_zone = factor(train_zone, levels = 1:4, labels = zone_labels),
    test_zone  = factor(test_zone,  levels = 1:4, labels = zone_labels),
    Method     = factor(Method, levels = c("LM", "OC4")),
    # make diagonal empty
    MSE        = ifelse(train_zone == test_zone, NA, MSE)
  )

#-----------------------------
# 3. Heatmap (faceted by method)
#-----------------------------
p <- ggplot(heat_long,
            aes(x = test_zone, y = train_zone, fill = MSE)) +
  geom_tile(color = "grey70") +
  scale_fill_gradient(
    low = "white",
    high = "#b2182b",
    na.value = "white",
    name = "MSE"
  ) +
  geom_text(
    aes(label = ifelse(is.na(MSE), "", sprintf("%.3f", MSE))),
    size = 4,
    color = "black",
    fontface = "bold"
  ) +
  labs(
    x = "Test Zone",
    y = "Train Zone",
    fill = "MSE"
  ) +
  coord_equal() +
  facet_wrap(~ Method) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 14),
    axis.title.x    = element_text(size = 16),
    axis.title.y    = element_text(size = 16),
    strip.text      = element_text(size = 16, face = "bold"),
    legend.title    = element_text(size = 15),
    legend.text     = element_text(size = 14)
  )

save_file <- "output/figure/plot_lat-abs-segment_heat-map_16_LM_OC4.pdf"
ggsave(save_file, plot = p, width = 8, height = 5)
