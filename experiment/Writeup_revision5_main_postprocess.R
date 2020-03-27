rm(list=ls())
load("../results/step5_clustering.RData")

esvd_curves$lineages

###################

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(100/255, 100/255, 100/255, alpha)) #gray
}

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})


num_order_vec <- c(5, rep(3,2), c(6,1,1,1,6,4), rep(2,2),  rep(5,2))
col_vec <- color_func(1)[num_order_vec]
col_vec3 <- color_func(0.3)[num_order_vec]
col_name <- c("orange", rep("bluish green", 2), c("bluish green", "yellow", "yellow", "yellow", "bluish green", "blue"), rep("skyblue", 2), rep("orange", 2))
order_vec <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                       order = order_vec,
                       col_name = col_name,
                       col_code = col_vec)
col_info

mean_vec <- t(sapply(1:13, function(i){
  idx <- which(cluster_labels == i)
  colMeans(esvd_embedding$u_mat[idx,])
}))

combn_mat <- combn(3,2)
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup5_main_2dplots_", k, ".png"),
      height = 1750, width = 1750, res = 300,
      units = "px")
  plot(x = esvd_embedding$u_mat[,i], y = esvd_embedding$u_mat[,j],
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       col = col_vec3[cluster_labels], pch = 16,
       main = "eSVD embedding and trajectories\n(Curved Gaussian)")


  curves <- esvd_curves$curves
  for(ll in 1:length(curves)){
    ord <- curves[[ll]]$ord
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 8)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black",
          lty = 3, lwd = 2)
  }

  text(mean_vec[,i], mean_vec[,j], labels = as.character(col_info$order), col = "black", cex = 1, font = 2)
  graphics.off()
}

#####################

# training vs testing

png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup5_main_esvd_training_testing.png"),
    height = 1200, width = 2000, res = 300,
    units = "px")
par(mfrow = c(1,2))
idx <- which.min(abs(esvd_angle_vec - 45))

nat_mat_list <- lapply(1:length(esvd_missing_list[[idx]]), function(i){
  esvd_missing_list[[idx]][[i]]$u_mat %*% t(esvd_missing_list[[idx]][[i]]$v_mat)
})

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 family = "curved_gaussian",
                                 scalar = paramMat[idx, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 max_points = 1e6)


plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 family = "curved_gaussian",
                                 scalar = paramMat[idx, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Testing set)")

graphics.off()

#### now for svd

load("../results/step2_naive_svd.RData")

png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup5_main_svd_training_testing.png"),
    height = 1200, width = 2000, res = 300,
    units = "px")
par(mfrow = c(1,2))
nat_mat_list <- lapply(1:length(svd_missing), function(i){
  svd_missing[[i]]$u %*% diag(svd_missing[[i]]$d) %*% t(svd_missing[[i]]$v)
})

tmp_mat <- do.call(rbind, lapply(1:length(nat_mat_list), function(i){
  cbind(dat_impute[missing_idx_list[[i]]], nat_mat_list[[i]][missing_idx_list[[i]]])
}))
sd_val <- sd(tmp_mat[,1] - tmp_mat[,2])

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 family = "gaussian", scalar = sd_val,
                                 main = "SVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 max_points = 1e6)


plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 family = "gaussian", scalar = sd_val,
                                 main = "SVD embedding:\nMatrix-completion diagnostic\n(Testing set)")

graphics.off()
