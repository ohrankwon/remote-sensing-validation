#===== 1-nearest-neighbour prediction in Rrs space
predict_1nn <- function(x_train, y_train, x0) {
  x_train <- as.matrix(x_train)
  x0 <- as.matrix(x0)

  m <- nrow(x0)
  y_pred <- numeric(m)

  for (i in 1:m) {
    d <- sqrt(rowSums((x_train - matrix(x0[i, ], nrow(x_train), ncol(x_train), byrow = TRUE))^2))
    nn_idx <- which.min(d)
    y_pred[i] <- y_train[nn_idx]
  }

  return(y_pred)
}

#===== spatiotemporal clustering (single linkage on the d_infinity dissimilarity)
spatiotemporal_cluster <- function(lat, lon, date, T_thresh = 7, S_thresh = 50) {
  if (length(lat) != length(lon) || length(lat) != length(date)) {
    stop("lat, lon, and date must have the same length.")
  }

  time <- as.Date(date)
  time_num <- as.numeric(time)

  delta_t <- abs(outer(time_num, time_num, "-"))

  coords <- cbind(lon, lat)
  delta_s <- geosphere::distm(coords, coords, fun = geosphere::distHaversine) / 1000

  d_inf <- pmax(delta_t / T_thresh, delta_s / S_thresh)
  d_inf_dist <- as.dist(d_inf)

  hc <- hclust(d_inf_dist, method = "single")
  clusters <- cutree(hc, h = 1)

  cluster_sizes <- as.vector(table(clusters))

  list(
    n_clusters           = length(unique(clusters)),
    cluster_size_summary = summary(cluster_sizes),
    cluster_size_table   = table(clusters),
    clusters             = clusters,
    hclust_object        = hc
  )
}
