rm(list=ls())

load("../results/step2_naive_svd.RData")

starting_lambda <- min(sapply(1:length(missing_idx_list), function(i){
  dat_NA <- dat_impute
  dat_NA[missing_idx_list[[i]]] <- NA

  softImpute::lambda0(dat_NA)
}))

lambda_vec <- seq(1, starting_lambda, length.out = 100)

quality_vec <- sapply(lambda_vec, function(lambda){
  print(lambda)

  k <- 5
  tmp_mat <- do.call(rbind, lapply(1:length(missing_idx_list), function(i){
    dat_NA <- dat_impute
    dat_NA[missing_idx_list[[i]]] <- NA

    softImpute_embedding <- softImpute::softImpute(dat_NA, rank.max = k, lambda = lambda)
    softImpute_pred <- softImpute_embedding$u %*% diag(softImpute_embedding$d) %*% t(softImpute_embedding$v)
    softImpute_pred <- softImpute_pred[missing_idx_list[[i]]]
    obs_val <- dat_impute[missing_idx_list[[i]]]

    cbind(obs_val, softImpute_pred)
  }))

  pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
  rad <- 2/5*max(tmp_mat[,1])
  ang <- as.numeric(acos(abs(c(0,1) %*% pca_res$rotation[,1])))
  ang * 180/pi
})

lambda <- lambda_vec[which.min(abs(quality_vec - 45))]

softImpute_missing_list <- lapply(1:length(missing_idx_list), function(i){
  dat_NA <- dat_impute
  dat_NA[missing_idx_list[[i]]] <- NA

  softImpute::softImpute(dat_NA, rank.max = k, lambda = lambda)
})

save(softImpute_missing_list, file = "../results/experiment.RData")

###############

load("../results/step5_clustering.RData")
load("../results/experiment.RData")


png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup5_main_svd_training_testing.png"),
    height = 1750, width = 3250, res = 300,
    units = "px")
par(mfrow = c(1,2))

nat_mat_list <- lapply(1:length(softImpute_missing_list), function(i){
  softImpute_missing_list[[i]]$u %*% diag(softImpute_missing_list[[i]]$d) %*% t(softImpute_missing_list[[i]]$v)
})

tmp_mat <- do.call(rbind, lapply(1:length(nat_mat_list), function(i){
  cbind(dat_impute[missing_idx_list[[i]]], nat_mat_list[[i]][missing_idx_list[[i]]])
}))
scalar <- sd(tmp_mat[,1] - tmp_mat[,2])

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 family = "gaussian", scalar = scalar,
                                 main = "SVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 max_points = 1e6)


plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 family = "gaussian", scalar = scalar,
                                 main = "SVD embedding:\nMatrix-completion diagnostic\n(Testing set)")

graphics.off()


