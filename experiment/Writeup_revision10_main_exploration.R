rm(list=ls())
load("../results/step4_factorization_cg_spca-vst_before_rescaling.RData")

cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(100/255, 100/255, 100/255, alpha)) #gray
}


num_order_vec_svd <- c(5, rep(3,2), c(1,1,1,1,1,1), rep(2,2),  rep(5,2))
col_vec_svd <- color_func(1)[num_order_vec_svd]
col_vec2_svd <- color_func(0.5)[num_order_vec_svd]
col_name_svd <- c("orange", rep("bluish green", 2), rep("yellow", 6), rep("skyblue", 2), rep("orange", 2))
order_vec_svd <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
col_info_svd <- data.frame(name = levels(cell_type_vec),
                           idx = sort(unique(cluster_labels)),
                           order = order_vec_svd,
                           col_name = col_name_svd,
                           col_code = col_vec_svd)
col_info_svd$factor_idx <- as.numeric(as.factor(col_info_svd$col_name))
col_info_svd[,c(5,6)] <- col_info_svd[,c(6,5)]
colnames(col_info_svd)[c(5,6)] <- colnames(col_info_svd)[c(6,5)]
col_info_svd
plotting_order_svd <- c(2,3,1,4)
cluster_center1 <- .compute_cluster_center(esvd_embedding$u_mat, .construct_cluster_matrix(cluster_labels))


combn_mat <- combn(3,2)

par(mfrow = c(1,3))
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]
  plot(x = esvd_embedding$u_mat[,i], y =  esvd_embedding$u_mat[,j],
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "SVD embedding and trajectories\n(Constant-variance Gaussian)",
       pch = 16, col = col_info_svd$col_code[cluster_labels])

  for(ll in 1:nrow(cluster_center1)){
    points(cluster_center1[ll,i], cluster_center1[ll,j], pch = 16, cex = 2.25, col = "black")
    points(cluster_center1[ll,i], cluster_center1[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }
}

################3
idx_choice <- 1
nat_mat_list <- lapply(1:length(esvd_missing_list[[idx_choice]]), function(i){
  esvd_missing_list[[idx_choice]][[i]]$u_mat %*% t(esvd_missing_list[[idx_choice]][[i]]$v_mat)
})

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

par(mfrow = c(1,2))
plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 family = fitting_distr,
                                 scalar = paramMat_esvd[idx_choice, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 max_points = 1e6)


plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 family = fitting_distr,
                                 scalar = paramMat_esvd[idx_choice, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Testing set)")

par(mfrow = c(1,2))
plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 family = fitting_distr,
                                 scalar = paramMat_esvd[idx_choice, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 max_points = 1e6, xlim = c(0, 100), ylim = c(0, 100))


plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 family = fitting_distr,
                                 scalar = paramMat_esvd[idx_choice, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Testing set)",
                                 xlim = c(0, 100), ylim = c(0, 100))

# training testing

dat_org <- log2(dat_impute/rescaling_factor+1)
nat_mat_list <- lapply(1:cv_trials, function(i){
  svd_missing_list[[i]]$u %*% diag(svd_missing_list[[i]]$d) %*% t(svd_missing_list[[i]]$v)
})
tmp_mat <- do.call(rbind, lapply(1:cv_trials, function(i){
  cbind(as.numeric(dat_org), as.numeric(nat_mat_list[[i]]))[missing_idx_list[[i]],]
}))
sd_val <- sd(tmp_mat[,1] - tmp_mat[,2])

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

par(mfrow = c(1,2))
plot_prediction_against_observed(dat_org, nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 family = "gaussian",
                                 scalar = sd_val,
                                 main = "SVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 max_points = 1e6)


plot_prediction_against_observed(dat_org, nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 family = "gaussian",
                                 scalar = sd_val,
                                 main = "SVD embedding:\nMatrix-completion diagnostic\n(Testing set)")

