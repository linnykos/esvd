var <- ls()
rm(list = var[var != "suffix"])
load(paste0("../results/step6_figures", suffix, ".RData"))

# make different variants of the 2D embedding plots
# 2D plots with curves (color)
# 2D plots in black and white
# 2D plots all combined in one plot (color)
# 2D plots without curves (color)

combn_mat <- combn(3,2)

######################

# esvd plots

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  png(filename = paste0("../../esvd_results/figure/main/esvd_2dplots_", k, ".png"),
      height = 1500, width = 1500, res = 300,
      units = "px")
  plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding and trajectories\n(Curved Gaussian)")

  for(ll in plotting_order_esvd) {
    target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx == ll)]
    idx <- which(cluster_labels %in% target_indices)
    points(x = esvd_embedding$u_mat[idx,i], y = esvd_embedding$u_mat[idx,j], pch = 16,
           col = col_info_esvd$col_code[target_indices[1]])
  }

  for(ll in 1:nrow(cluster_center_esvd)){
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
  }


  curves <- esvd_curves$curves
  for(ll in 1:length(curves)) {
    ord <- curves[[ll]]$ord
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 15)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = col_vec_short[ll], lwd = 7)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black",
          lty = 3, lwd = 3)
  }

  graphics.off()
}

png(filename = paste0("../../esvd_results/figure/main/esvd_2dplots.png"),
    height = 830, width = 2300, res = 300,
    units = "px")
par(mfrow = c(1,3), mar = c(4,4,4,1))
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]
  plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = ifelse(k==2, "eSVD embedding and trajectories\n(Curved Gaussian)", ""))

  for(ll in plotting_order_esvd) {
    target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx == ll)]
    idx <- which(cluster_labels %in% target_indices)
    points(x = esvd_embedding$u_mat[idx,i], y = esvd_embedding$u_mat[idx,j], pch = 16,
           col = col_vec2_esvd[target_indices[1]])
  }

  for(ll in 1:nrow(cluster_center_esvd)){
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2, col = "black")
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
  }


  curves <- esvd_curves$curves
  for(ll in 1:length(curves)) {
    ord <- curves[[ll]]$ord
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 8)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = col_vec_short[ll], lwd = 5)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black",
          lty = 3, lwd = 2)
  }
}
graphics.off()

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  png(filename = paste0("../../esvd_results/figure/main/esvd_2dplots_", k, "_nocurve.png"),
      height = 1500, width = 1500, res = 300,
      units = "px")
  plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding and trajectories\n(Curved Gaussian)")

  for(ll in plotting_order_esvd) {
    target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx == ll)]
    idx <- which(cluster_labels %in% target_indices)
    points(x = esvd_embedding$u_mat[idx,i], y = esvd_embedding$u_mat[idx,j], pch = 16,
           col = col_info_esvd$col_code[target_indices[1]])
  }

  for(ll in 1:nrow(cluster_center_esvd)){
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
  }

  graphics.off()
}

########################################

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  png(filename = paste0("../../esvd_results/figure/main/svd_2dplots_", k, ".png"),
      height = 1500, width = 1500, res = 300,
      units = "px")
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

  graphics.off()
}

png(filename = paste0("../../esvd_results/figure/main/svd_2dplots.png"),
    height = 830, width = 2300, res = 300,
    units = "px")
par(mfrow = c(1,3), mar = c(4,4,4,1))

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  plot(NA, xlim = range(svd_embedding[,i]), ylim = range(svd_embedding[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = ifelse(k == 2, "SVD embedding and trajectories\n(Constant-variance Gaussian)", ""))

  for(ll in plotting_order_svd) {
    target_indices <- col_info_svd$idx[which(col_info_svd$factor_idx == ll)]
    idx <- which(cluster_labels %in% target_indices)
    points(x = svd_embedding[idx,i], y = svd_embedding[idx,j], pch = 16,
           col = col_vec2_svd[target_indices[1]])
  }


  for(ll in 1:nrow(cluster_center_svd)){
    points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 2, col = "black")
    points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }


  curves <- svd_curves$curves
  for(ll in 1:length(curves)) {
    ord <- curves[[ll]]$ord
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 5)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black", lwd = 3)
  }
}
graphics.off()


for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  png(filename = paste0("../../esvd_results/figure/main/svd_2dplots_", k, "_nocurve.png"),
      height = 1500, width = 1500, res = 300,
      units = "px")
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

  graphics.off()
}
