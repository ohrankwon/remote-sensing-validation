#===== workspace
# Groups the matched observations into spatiotemporal clusters using single
# linkage on the d_infinity dissimilarity with T = 7 days and S = 50 km
# (Lemma 1 of the main paper), cutting the dendrogram at height 1.
#   in:  data/chl_data_log_source.csv
#   out: data/chl_data_log_cluster.csv
rm(list = ls())

library(geosphere)
source("R/functions.R")

#===== data
data <- read.csv("data/chl_data_log_source.csv", header = T)

result <- spatiotemporal_cluster(
  lat = data$lat,
  lon = data$lon,
  date = data$date,
  T_thresh = 7,
  S_thresh = 50
)

# combine clusters info on data
data$clusters <- result$clusters
write.csv(data,
          file = "data/chl_data_log_cluster.csv",
          row.names = FALSE)

# Total number of clusters (8,177 as reported in Section 2.3 of the main paper)
result$n_clusters

# Summary of cluster sizes
result$cluster_size_summary
