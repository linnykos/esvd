var <- ls()
rm(list = var[var != "suffix"])
load(paste0("../results/step7_figures", suffix, ".RData"))

# make different variants of the 2D embedding plots
# 2D plots with curves (color)
# 2D plots in black and white
# 2D plots all combined in one plot (color)
# 2D plots without curves (color)

######################

# esvd plots

esvd_curves_prepared <- lapply(esvd_curves_short$curves, function(curve){
  prepare_trajectory(curve, target_length = 11)
})

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  grDevices::png(filename = paste0("../../esvd_results/figure/main/esvd_2dplots_", k, ".png"),
      height = 1500, width = 1500, res = 300,
      units = "px")
  graphics::plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding and trajectories\n(Curved Gaussian)")

  for(ll in plotting_order_esvd) {
    target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx %in% ll)]
    idx <- which(cluster_labels %in% target_indices)
    graphics::points(x = esvd_embedding$u_mat[idx,i], y = esvd_embedding$u_mat[idx,j], pch = 16,
           col = col_vec2_esvd[cluster_labels[idx]])
  }

  curves <- esvd_curves_prepared
  for(ll in rev(1:length(curves))) {
    graphics::lines(x = curves[[ll]][, i], y = curves[[ll]][, j], col = "white", lwd = 15)
    graphics::lines(x = curves[[ll]][, i], y = curves[[ll]][, j], col = col_vec_short[ll], lwd = 7)
    graphics::lines(x = curves[[ll]][, i], y = curves[[ll]][, j], col = "black",
          lty = 3, lwd = 3)
  }

  for(ll in 1:nrow(cluster_center_esvd)){
    graphics::points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    graphics::points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
  }

  grDevices::graphics.off()
}

grDevices::png(filename = paste0("../../esvd_results/figure/main/esvd_2dplots.png"),
    height = 830, width = 2300, res = 300,
    units = "px")
graphics::par(mfrow = c(1,3), mar = c(4,4,4,1))
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]
  graphics::plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = ifelse(k==2, "eSVD embedding and trajectories\n(Curved Gaussian)", ""))

  for(ll in plotting_order_esvd) {
    target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx %in% ll)]
    idx <- which(cluster_labels %in% target_indices)
    graphics::points(x = esvd_embedding$u_mat[idx,i], y = esvd_embedding$u_mat[idx,j], pch = 16,
           col = col_vec2_esvd[cluster_labels[idx]])
  }

  curves <- esvd_curves_prepared
  for(ll in rev(1:length(curves))) {
    graphics::lines(x = curves[[ll]][, i], y = curves[[ll]][, j], col = "white", lwd = 15)
    graphics::lines(x = curves[[ll]][, i], y = curves[[ll]][, j], col = col_vec_short[ll], lwd = 7)
    graphics::lines(x = curves[[ll]][, i], y = curves[[ll]][, j], col = "black",
          lty = 3, lwd = 3)
  }

  for(ll in 1:nrow(cluster_center_esvd)){
    graphics::points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    graphics::points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
  }
}
grDevices::graphics.off()

grDevices::png(filename = paste0("../../esvd_results/figure/main/esvd_2dplots_nocurve.png"),
               height = 830, width = 2300, res = 300,
               units = "px")
graphics::par(mfrow = c(1,3), mar = c(4,4,4,1))
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]
  graphics::plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
                 asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
                 main = ifelse(k==2, "eSVD embedding and trajectories\n(Curved Gaussian)", ""))

  for(ll in plotting_order_esvd) {
    target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx %in% ll)]
    idx <- which(cluster_labels %in% target_indices)
    graphics::points(x = esvd_embedding$u_mat[idx,i], y = esvd_embedding$u_mat[idx,j], pch = 16,
                     col = col_vec2_esvd[cluster_labels[idx]])
  }

  for(ll in 1:nrow(cluster_center_esvd)){
    graphics::points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    graphics::points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
  }
}
grDevices::graphics.off()

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  grDevices::png(filename = paste0("../../esvd_results/figure/main/esvd_2dplots_", k, "_nocurve.png"),
      height = 1500, width = 1500, res = 300,
      units = "px")
  graphics::plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding and trajectories\n(Curved Gaussian)")

  for(ll in plotting_order_esvd) {
    target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx %in% ll)]
    idx <- which(cluster_labels %in% target_indices)
    graphics::points(x = esvd_embedding$u_mat[idx,i], y = esvd_embedding$u_mat[idx,j], pch = 16,
           col = col_vec2_esvd[cluster_labels[idx]])
  }

  for(ll in 1:nrow(cluster_center_esvd)){
    graphics::points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    graphics::points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
  }

  grDevices::graphics.off()
}

########################################

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  grDevices::png(filename = paste0("../../esvd_results/figure/main/svd_2dplots_", k, ".png"),
      height = 1500, width = 1500, res = 300,
      units = "px")
  graphics::plot(NA, xlim = range(svd_embedding[,i]), ylim = range(svd_embedding[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "SVD embedding and trajectories\n(Constant-variance Gaussian)")

  for(ll in plotting_order_svd) {
    target_indices <- col_info_svd$idx[which(col_info_svd$factor_idx %in% ll)]
    idx <- which(cluster_labels %in% target_indices)
    graphics::points(x = svd_embedding[idx,i], y = svd_embedding[idx,j], pch = 16,
           col = col_vec2_svd[cluster_labels[idx]])
  }

  curves <- svd_curves_short$curves
  for(ll in 1:length(curves)) {
    ord <- curves[[ll]]$ord
    graphics::lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 8)
    graphics::lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black", lwd = 5)
  }

  for(ll in 1:nrow(cluster_center_svd)){
    graphics::points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 2, col = "black")
    graphics::points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }


  grDevices::graphics.off()
}

grDevices::png(filename = paste0("../../esvd_results/figure/main/svd_2dplots.png"),
    height = 830, width = 2300, res = 300,
    units = "px")
graphics::par(mfrow = c(1,3), mar = c(4,4,4,1))

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  graphics::plot(NA, xlim = range(svd_embedding[,i]), ylim = range(svd_embedding[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = ifelse(k == 2, "SVD embedding and trajectories\n(Constant-variance Gaussian)", ""))

  for(ll in plotting_order_svd) {
    target_indices <- col_info_svd$idx[which(col_info_svd$factor_idx %in% ll)]
    idx <- which(cluster_labels %in% target_indices)
    graphics::points(x = svd_embedding[idx,i], y = svd_embedding[idx,j], pch = 16,
           col = col_vec2_svd[cluster_labels[idx]])
  }

  curves <- svd_curves_short$curves
  for(ll in 1:length(curves)) {
    ord <- curves[[ll]]$ord
    graphics::lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 5)
    graphics::lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black", lwd = 3)
  }

  for(ll in 1:nrow(cluster_center_svd)){
    graphics::points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 2, col = "black")
    graphics::points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }
}
grDevices::graphics.off()



grDevices::png(filename = paste0("../../esvd_results/figure/main/svd_2dplots_nocurve.png"),
               height = 830, width = 2300, res = 300,
               units = "px")
graphics::par(mfrow = c(1,3), mar = c(4,4,4,1))

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  graphics::plot(NA, xlim = range(svd_embedding[,i]), ylim = range(svd_embedding[,j]),
                 asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
                 main = ifelse(k == 2, "SVD embedding and trajectories\n(Constant-variance Gaussian)", ""))

  for(ll in plotting_order_svd) {
    target_indices <- col_info_svd$idx[which(col_info_svd$factor_idx %in% ll)]
    idx <- which(cluster_labels %in% target_indices)
    graphics::points(x = svd_embedding[idx,i], y = svd_embedding[idx,j], pch = 16,
                     col = col_vec2_svd[cluster_labels[idx]])
  }

  for(ll in 1:nrow(cluster_center_svd)){
    graphics::points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 2, col = "black")
    graphics::points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }
}
grDevices::graphics.off()

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  grDevices::png(filename = paste0("../../esvd_results/figure/main/svd_2dplots_", k, "_nocurve.png"),
      height = 1500, width = 1500, res = 300,
      units = "px")
  graphics::plot(NA, xlim = range(svd_embedding[,i]), ylim = range(svd_embedding[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "SVD embedding and trajectories\n(Constant-variance Gaussian)")

  for(ll in plotting_order_svd) {
    target_indices <- col_info_svd$idx[which(col_info_svd$factor_idx %in% ll)]
    idx <- which(cluster_labels %in% target_indices)
    graphics::points(x = svd_embedding[idx,i], y = svd_embedding[idx,j], pch = 16,
           col = col_vec2_svd[cluster_labels[idx]])
  }

  for(ll in 1:nrow(cluster_center_svd)){
    graphics::points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 2, col = "black")
    graphics::points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }

  grDevices::graphics.off()
}

