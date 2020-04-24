par(mfrow = c(1,2))
nat_mat_list <- lapply(1:length(svd_missing), function(i){
  svd_missing[[i]]$u %*% diag(svd_missing[[i]]$d) %*% t(svd_missing[[i]]$v)
})

tmp_mat <- do.call(rbind, lapply(1:length(nat_mat_list), function(i){
  cbind(log2(dat/rescaling_factor + 1)[missing_idx_list[[i]]], nat_mat_list[[i]][missing_idx_list[[i]]])
}))
sd_val <- sd(tmp_mat[,1] - tmp_mat[,2])

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

plot_prediction_against_observed(log2(dat_impute+1), nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 family = "gaussian", scalar = sd_val,
                                 main = "SVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 max_points = 1e6)


plot_prediction_against_observed(log2(dat_impute+1), nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 family = "gaussian", scalar = sd_val,
                                 main = "SVD embedding:\nMatrix-completion diagnostic\n(Testing set)")

############################

par(mfrow = c(1,2))

nat_mat_list <- lapply(1:length(esvd_missing_list[[esvd_angle_res$idx]]), function(i){
  esvd_missing_list[[esvd_angle_res$idx]][[i]]$u_mat %*% t(esvd_missing_list[[esvd_angle_res$idx]][[i]]$v_mat)
})

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 family = "neg_binom",
                                 scalar = paramMat_esvd[esvd_angle_res$idx, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 max_points = 1e6)


plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 family = "neg_binom",
                                 scalar = paramMat_esvd[esvd_angle_res$idx, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Testing set)")
