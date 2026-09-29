#===== workspace
rm(list = ls())

library(matrixStats) # rowMaxs
library(dplyr)
source("R/functions.R")

#===== set-up
# All 16 train-test combinations of the 4 latitude zones. Row order follows
# expand.grid(train_zone, test_zone).
num_repeat <- 16
set_seed <- 513971

#===== data
data = data_org <- read.csv("data/chl_data_log_cluster.csv", header = T)
names(data)[4:10] = names(data_org)[4:10] <- c(paste0("x",1:6), "y")

#===== functions
mse_se <- function(yhat, y){
  mse = sum( (yhat - y)^2 ) / length(y)
  se = sd( (yhat - y)^2 ) / sqrt(length(y))

  return( c(mse, se) )
}

mse_se_diff <- function(yhat, y, yhat_rel){
  mse = sum( (yhat - y)^2 - (yhat_rel - y)^2 ) / length(y)
  se = sd( (yhat - y)^2 - (yhat_rel - y)^2 ) / sqrt(length(y))

  return( c(mse, se) )
}

cut_off <- c(0, 23.4, 35, 60, 90)  # 4 latitude ranges
n_zone  <- 4

# all train-test combinations: (1,1), (1,2), ..., (4,4)  -> 16 rows
pairs <- expand.grid(train_zone = 1:n_zone,
                     test_zone  = 1:n_zone)

#===== analysis
for(i in 1:num_repeat) {

  #===== output format
  if (i==1) {
    set.seed(set_seed)
    eval = eval_diff <- matrix(NA, num_repeat, 16)
    colnames(eval) = colnames(eval_diff) <- c("nn_mse", "lm_mse", "oc4_mse",
                                              "nn_se", "lm_se", "oc4_se",
                                              "n_train", "n_test",
                                              "train_lat_l", "train_lat_u", "test_lat_l", "test_lat_u",
                                              "n_sub_test", "oc4_sub_se", "nn_sub_se", "lm_sub_se")
  }

  #===== train/test data
  data <- data_org
  unique_cluster <- unique(data$clusters)

  cluster_lat_mean <- tapply(data$lat, data$clusters, mean)
  cluster_lat_mean <- abs(cluster_lat_mean)

  # split into train/test
  tr_z <- pairs$train_zone[i]
  te_z <- pairs$test_zone[i]
  train_clusters <- which(cluster_lat_mean >= cut_off[tr_z] &
                            cluster_lat_mean <  cut_off[tr_z + 1])
  test_clusters  <- which(cluster_lat_mean >= cut_off[te_z] &
                            cluster_lat_mean <  cut_off[te_z + 1])

  ## indices in the data
  train_id <- which(data$clusters %in% train_clusters)
  test_id  <- which(data$clusters %in% test_clusters)

  ## X-type variables
  train    <- data[train_id, 4:10]
  test     <- data[test_id,  4:10]

  ## ST-type variables
  test_st_cluster <- data[ test_id , c(1:3,10,12)]

  test_st_cluster_with_id <- test_st_cluster %>%
    mutate(row_id = row_number())  # Keep original row index

  sample_per_cluster <- test_st_cluster_with_id %>%
    group_by(clusters) %>%
    slice_sample(n = 1) %>%
    ungroup()

  kk <- sample_per_cluster$row_id

  eval[i, c("n_train","n_test")] = eval_diff[i, c("n_train","n_test")] <- c(dim(train)[1], dim(test)[1])
  eval[i, c("n_sub_test")] = eval_diff[i, c("n_sub_test")] <- length(kk)

  # latitude bounds for each TRAIN range
  train_lat_l <- min(cluster_lat_mean[train_clusters])
  train_lat_u <- max(cluster_lat_mean[train_clusters])

  eval[i, c("train_lat_l","train_lat_u")] <- c(train_lat_l, train_lat_u)
  eval_diff[i, c("train_lat_l","train_lat_u")] <- c(train_lat_l, train_lat_u)

  # latitude bounds for each TEST range
  test_lat_l <- min(cluster_lat_mean[test_clusters])
  test_lat_u <- max(cluster_lat_mean[test_clusters])

  eval[i, c("test_lat_l","test_lat_u")] <- c(test_lat_l, test_lat_u)
  eval_diff[i, c("test_lat_l","test_lat_u")] <- c(test_lat_l, test_lat_u)

  # OC4
  r <- rowMaxs( as.matrix(train[,2:4]) ) - train[,5]
  train_r <- cbind(r, r^2, r^3, r^4)
  r <- rowMaxs( as.matrix(test[,2:4]) ) - test[,5]
  test_r <- cbind(r, r^2, r^3, r^4)
  fit_oc4 <- lm(train[, 7] ~ ., data = as.data.frame(train_r))
  yhat_oc4 <- predict(fit_oc4, newdata = as.data.frame(test_r))
  eval_oc4 <- mse_se(yhat = yhat_oc4, y = test[,7])
  eval[i, c("oc4_mse","oc4_se")] <- eval_oc4
  eval_oc4 <- mse_se_diff(yhat = yhat_oc4, y = test[,7], yhat_rel = yhat_oc4)
  eval_diff[i, c("oc4_mse","oc4_se")] <- eval_oc4

  eval_oc4 <- mse_se(yhat = yhat_oc4[kk], y = (test[,7])[kk] )
  eval[i, c("oc4_sub_se")] <- eval_oc4[2]
  eval_oc4 <- mse_se_diff(yhat = yhat_oc4[kk], y = (test[,7])[kk], yhat_rel = yhat_oc4[kk])
  eval_diff[i, c("oc4_sub_se")] <- eval_oc4[2]

  # 1NN
  yhat_sup <- predict_1nn(x_train = train[,-7], y_train = train[,7], test[,-7])
  eval_sup <- mse_se(yhat = yhat_sup, y = test[,7])
  eval[i, c("nn_mse","nn_se")] <- eval_sup
  eval_sup <- mse_se_diff(yhat = yhat_sup, y = test[,7], yhat_rel = yhat_oc4)
  eval_diff[i, c("nn_mse","nn_se")] <- eval_sup

  eval_sup <- mse_se(yhat = yhat_sup[kk], y = (test[,7])[kk] )
  eval[i, c("nn_sub_se")] <- eval_sup[2]
  eval_sup <- mse_se_diff(yhat = yhat_sup[kk], y = (test[,7])[kk], yhat_rel = yhat_oc4[kk])
  eval_diff[i, c("nn_sub_se")] <- eval_sup[2]

  # LM
  fit_lm <- lm(train[, 7] ~ ., data = as.data.frame(train[,-7]))
  yhat_lm <- predict(fit_lm, newdata = as.data.frame(test[,-7]))
  eval_lm <- mse_se(yhat = yhat_lm, y = test[,7])
  eval[i, c("lm_mse","lm_se")] <- eval_lm
  eval_lm <- mse_se_diff(yhat = yhat_lm, y = test[,7], yhat_rel = yhat_oc4)
  eval_diff[i, c("lm_mse","lm_se")] <- eval_lm

  eval_lm <- mse_se(yhat = yhat_lm[kk], y = (test[,7])[kk] )
  eval[i, c("lm_sub_se")] <- eval_lm[2]
  eval_lm <- mse_se_diff(yhat = yhat_lm[kk], y = (test[,7])[kk], yhat_rel = yhat_oc4[kk])
  eval_diff[i, c("lm_sub_se")] <- eval_lm[2]

}

file.name <- "output/RData/analysis_cross_region_grid.RData"
save(eval, eval_diff, file = file.name)
