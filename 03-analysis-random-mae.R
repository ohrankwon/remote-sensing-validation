#===== workspace
rm(list = ls())

library(matrixStats) # rowMaxs
library(dplyr)
library(withr)
source("R/functions.R")

#===== set-up
# Identical to 02-analysis-random-mse.R except that the error metric is the mean
# ABSOLUTE error rather than the mean squared error (see mse_se / mse_se_diff
# below). Used for the robustness check in Figure S4.
# The script runs BOTH splitting strategies and writes one RData file for each:
#   st_cluster = TRUE  -> "Clust":   whole clusters are assigned to train/test
#   st_cluster = FALSE -> "NoClust": every observation is made its own cluster,
#                                    so the split degenerates to the standard
#                                    observation-level random split.
cut_off_random <- 0.5 # proportion of samples used for training
num_repeat <- 1000
set_seed <- 513971

#===== data
data_full <- read.csv("data/chl_data_log_cluster.csv", header = T)
names(data_full)[4:10] <- c(paste0("x",1:6), "y")

#===== functions
# mean ABSOLUTE error and its standard error
mse_se <- function(yhat, y){
  mse = sum( abs(yhat - y) ) / length(y)
  se = sd( abs(yhat - y) ) / sqrt(length(y))

  return( c(mse, se) )
}

mse_se_diff <- function(yhat, y, yhat_rel){
  mse = sum( abs(yhat - y) - abs(yhat_rel - y) ) / length(y)
  se = sd( abs(yhat - y) - abs(yhat_rel - y) ) / sqrt(length(y))

  return( c(mse, se) )
}

#===== analysis
for (st_cluster in c(TRUE, FALSE)) {

  data = data_org <- data_full

  if (st_cluster == FALSE) {
    # every observation becomes its own cluster -> observation-level random split
    data$clusters = data_org$clusters <- 1:length(data$clusters)
  }

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
    # shuffle data
    data <- data_org
    unique_cluster <- unique(data$clusters)
    shuffled_cluster <- sample(unique_cluster)
    label_mapping <- setNames(shuffled_cluster, unique_cluster)
    data$cluster <- label_mapping[as.character(data$clusters)]
    data <- data[order(data$clusters),]

    cluster_date_loc_count <- table( as.factor(data$cluster) )
    cluster_date_loc_cumprop <- cumsum( cluster_date_loc_count ) / sum( cluster_date_loc_count )

    # split into train/test
    cut_off_cluster <- which(cluster_date_loc_cumprop >= cut_off_random)[1]

    train <- data[ which(data$cluster <= cut_off_cluster) , 4:10]
    test <- data[ which(data$cluster > cut_off_cluster) , 4:10]
    test_st_cluster <- data[ which(data$cluster > cut_off_cluster) , c(1:3,10,12,13)]

    test_st_cluster_with_id <- test_st_cluster %>%
      mutate(row_id = row_number())  # Keep original row index

    # drawn from an isolated RNG stream so it does not depend on the methods fitted
    sample_per_cluster <- with_seed(i*(set_seed-12345), {
      test_st_cluster_with_id %>%
        group_by(clusters) %>%
        slice_sample(n = 1) %>%
        ungroup()
    })

    kk <- sample_per_cluster$row_id

    eval[i, c("n_train","n_test")] = eval_diff[i, c("n_train","n_test")] <- c(dim(train)[1], dim(test)[1])
    eval[i, c("n_sub_test")] = eval_diff[i, c("n_sub_test")] <- length(kk)

    # OC4
    # variables: 'Rrs_412', 'Rrs_443', 'Rrs_490', 'Rrs_510', 'Rrs_560', 'Rrs_665'
    # R = log10( max( Rrs_443/Rrs_560, Rrs_490/Rrs_560, Rrs_510/Rrs_560 ) ). In fact, originally, max {Rrs 443; 490; 510} / Rrs 555
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

  file.name <- "output/RData/analysis_random_mae_clust.RData"
  if(st_cluster==FALSE){
    file.name <- "output/RData/analysis_random_mae_noclust.RData"
  }
  save(eval, eval_diff, file = file.name)

}
