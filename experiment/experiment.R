rm(list=ls())

load("../results/lingxue_analysis.RData")

sapply(fit_all_list, length)
sapply(fit_list, length)
sapply(fit_list[[1]]$neg_binom_missing, length)

dat_impute <- preprocessing_list[[1]]$dat_impute

sapply(1:7, function(i){
  nat_mat <- fit_list[[1]]$neg_binom_missing[[i]]$u_mat %*% t(fit_list[[1]]$neg_binom_missing[[i]]$v_mat)
  plot_prediction_against_observed(dat_impute, nat_mat, family = "neg_binom", missing_idx = missing_idx,
                                   scalar = neg_binom_vec[i], plot = F)
})

plot(fit_list[[1]]$neg_binom_missing[[i]]$u_mat[,1], fit_list[[1]]$neg_binom_missing[[i]]$u_mat[,2],
     asp = T, pch = 16, col = as.numeric(preprocessing_list[[1]]$label_vec))

sapply(1:7, function(i){
  nat_mat <- fit_list[[1]]$curved_gaussian_missing[[i]]$u_mat %*% t(fit_list[[1]]$curved_gaussian_missing[[i]]$v_mat)
  plot_prediction_against_observed(dat_impute, nat_mat, family = "curved_gaussian", missing_idx = missing_idx,
                                   scalar = curved_gaussian_vec[i], plot = F)
})
