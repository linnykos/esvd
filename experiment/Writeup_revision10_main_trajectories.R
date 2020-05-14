rm(list=ls())
load("../../esvd/results/step5_trajectory_original.RData")

p <- 3
set.seed(10)
esvd_curves <- eSVD::slingshot(esvd_embedding$u_mat[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                               cluster_group_list = cluster_group_list, shrink = 2,
                               verbose = T, upscale_factor = 1,
                               squared = T)


set.seed(10)
svd_curves <- slingshot(svd_embedding[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                        cluster_group_list = cluster_group_list,
                        verbose = T, upscale_factor = 1, shrink = 2, max_iter = 2,
                        squared = F)

#############################3

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(100/255, 100/255, 100/255, alpha)) #gray
}

combn_mat <- combn(3,2)

#####################

# info for the upcoming svd plots

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

cluster_center_svd <- .compute_cluster_center(svd_embedding[,1:3], .construct_cluster_matrix(cluster_labels))


# info for the upcoming esvd plots

num_order_vec_esvd <- c(5, rep(3,2), c(1,1,1,1,6,4), rep(2,2),  rep(5,2))
col_vec_esvd <- color_func(1)[num_order_vec_esvd]
col_vec2_esvd <- color_func(0.5)[num_order_vec_esvd]
col_name_esvd <- c("orange", rep("bluish green", 2), c("yellow", "yellow", "yellow", "yellow", "gray", "blue"), rep("skyblue", 2), rep("orange", 2))
order_vec_esvd <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
col_info_esvd <- data.frame(name = levels(cell_type_vec),
                            idx = sort(unique(cluster_labels)),
                            order = order_vec_esvd,
                            col_name = col_name_esvd,
                            col_code = col_vec_esvd)
col_info_esvd$factor_idx <- as.numeric(as.factor(col_info_esvd$col_name))
col_info_esvd[,c(5,6)] <- col_info_esvd[,c(6,5)]
colnames(col_info_esvd)[c(5,6)] <- colnames(col_info_esvd)[c(6,5)]
col_info_esvd
col_vec_short <- color_func(0.9)[c(1,4)]
plotting_order_esvd <- c(3,4,5,2,6,1)

cluster_center_esvd <- .compute_cluster_center(esvd_embedding$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))

###############################################

par(mfrow = c(1,3))

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  plot(NA, xlim = range(svd_embedding[,i]), ylim = range(svd_embedding[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "SVD embedding and trajectories\n(Constant-variance Gaussian)")

  for(ll in plotting_order_svd) {
    target_indices <- col_info_svd$idx[which(col_info_svd$factor_idx == ll)]
    idx <- which(cluster_labels %in% target_indices)
    points(x = svd_embedding[idx,i], y = svd_embedding[idx,j], pch = 16,
           col = col_info_svd$col_code[target_indices[1]])
  }


  for(ll in 1:nrow(cluster_center_svd)){
    points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 2, col = "black")
    points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }


  curves <- svd_curves$curves
  for(ll in 1:length(curves)) {
    ord <- curves[[ll]]$ord
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 8)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black", lwd = 5)
  }
}
