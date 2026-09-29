#===== workspace
rm(list = ls())

library(matrixStats) # rowMaxs
library(dplyr)
source("R/functions.R")

#===== set-up
# Combined cross-region + temporal scenarios, split at 2010-01-01.
#   M1 (i = 1): train on North Atlantic Mid-Latitude clusters before 2010,
#               test on ALL global clusters after 2010 (region + time shift).
#   M2 (i = 2): train on North Atlantic Mid-Latitude clusters before 2010,
#               test on North Atlantic Mid-Latitude clusters after 2010 (time shift only).
num_repeat <- 2
set_seed <- 513971

# North Atlantic Mid-Latitude box used for the training region
box_lat_l  <- 35
box_lat_u  <- 60
box_lon_l  <- -80
box_lon_u  <- 0
split_date <- as.numeric(as.Date("2010-01-01"))

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

#===== analysis
for(i in 1:num_repeat) {

  #===== output format
  if (i==1) {
    set.seed(set_seed)
    eval = eval_diff <- matrix(NA, num_repeat, 12)
    colnames(eval) = colnames(eval_diff) <- c("nn_mse", "lm_mse", "oc4_mse",
                                              "nn_se", "lm_se", "oc4_se",
                                              "n_train", "n_test",
                                              "n_sub_test", "oc4_sub_se", "nn_sub_se", "lm_sub_se")
  }

  #===== train/test data
  data <- data_org

  cluster_lat_mean  <- tapply(data$lat, data$clusters, mean)
  cluster_lon_mean  <- tapply(data$lon, data$clusters, mean)
  cluster_date_mean <- tapply(as.numeric(as.Date(data$date)), data$clusters, mean)

  in_box <- cluster_lat_mean >= box_lat_l & cluster_lat_mean <  box_lat_u &
            cluster_lon_mean >= box_lon_l & cluster_lon_mean <  box_lon_u

  # both scenarios train on the same clusters: in-box and before the split date
  train_clusters <- which(in_box & cluster_date_mean < split_date)

  if(i==1){
    # M1: test on all clusters after the split date
    test_clusters <- which(cluster_date_mean > split_date)
  }

  if(i==2){
    # M2: test only on in-box clusters after the split date
    test_clusters <- which(in_box & cluster_date_mean > split_date)
  }

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

file.name <- "output/RData/analysis_combined.RData"
save(eval, eval_diff, file = file.name)
