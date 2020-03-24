rm(list=ls())
load("../experiment/experiment.RData")

for(i in 1:nrow(paramMat)){
  nat_mat_list <- lapply(1:cv_trials, function(j){
    fit_all_list[[i]][[j]]$u_mat %*% t(fit_all_list[[i]][[j]]$v_mat)
  })

  eSVD::plot_prediction_against_observed(dat, nat_mat_list = nat_mat_list,
                                         scalar = paramMat[i, "scalar"],
                                         family = "neg_binom", missing_idx_list = missing_idx_list,
                                         plot = T, main = i)
}


i=1
j=1
cluster_labels <- rep(1:4, each = vec["n_each"])
plot(fit_all_list[[i]][[j]]$u_mat[,1], fit_all_list[[i]][[j]]$u_mat[,2], pch = 16, col = cluster_labels, asp = T)


####

i <- 4
cluster_labels <- rep(1:4, each = vec["n_each"])
plot(fit_list[[i]]$u_mat[,1], fit_list[[i]]$u_mat[,2], pch = 16, col = cluster_labels, asp = T)


missing_idx_list <- list(missing_idx)

quality_vec <- sapply(1:nrow(paramMat), function(i){
  nat_mat <-  fit_list[[i]]$u_mat %*% t(fit_list[[i]]$v_mat)
  nat_mat_list <- list(nat_mat)

  family = "neg_binom"
  scalar = paramMat[i, "scalar"]

  pred_mat_list <- lapply(nat_mat_list, function(nat_mat){
    compute_mean(nat_mat, family = family, scalar = scalar)
  })

  tmp_mat <- do.call(rbind, lapply(1:length(nat_mat_list), function(i){
    cbind(dat[missing_idx_list[[i]]], pred_mat_list[[i]][missing_idx_list[[i]]])
  }))

  idx <- which(tmp_mat[,1] != 0)
  tmp_mat <- tmp_mat[idx,]

  pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
  rad <- 2/5*max(tmp_mat[,1])
  ang <- as.numeric(acos(abs(c(0,1) %*% pca_res$rotation[,1])))

  ang * 180/pi
})
