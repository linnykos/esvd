var <- ls()
rm(list = var[var != "suffix"])
load(paste0("../results/step6_figures", suffix, ".RData"))


png(paste0("../../esvd_results/figure/main/svd_table.png"),
    height = 1300, width = 1500, res = 300, units = "px")
plot(NA,  asp = T, xlab = paste0("Latent dimension ", dim1), ylab = paste0("Latent dimension ", dim2),
     main = "SVD embedding\n(Constant-variance Gaussian)",
     xlim = c(-10,10), ylim = range(-y), axes = F)
points(x = x, y = -y, col = col_vec2_svd[cluster_labels], pch = 16)
axis(1); axis(2)

for(i in 3:6){
  if(i == 3){
    clust_idx <- which(col_info_svd[,"level"] %in% c(1:3))
  } else {
    clust_idx <- which(col_info_svd[,"level"] %in% i)
  }
  idx <- which(cluster_labels %in% clust_idx)
  kde_est <- MASS::kde2d(svd_embedding[idx,dim1], -svd_embedding[idx,dim2], n = 500, lims = c(range(x), range(-y)))
  col_val <- col_info[clust_idx[1], "col_code"]
  quant_level <- quantile(kde_est$z, probs = 0.925)
  contour(kde_est, add = T, drawlabels = F, lwd = 5, levels = quant_level, col = "white")
  contour(kde_est, add = T, drawlabels = F, lwd = 3, levels = quant_level, col = col_val)
  contour(kde_est, add = T, drawlabels = F, lwd = 1, levels = quant_level,  lty = 3, col = "black")
}
graphics.off()


###########################

# esvd density

png(paste0("../../esvd_results/figure/main/esvd_density.png"),
    height = 1500, width = 1500, res = 300, units = "px")
x <- esvd_embedding$u_mat[,dim1]
y <- esvd_embedding$u_mat[,dim2]
plot(NA,  asp = T, xlab = paste0("Latent dimension ", dim1), ylab = paste0("Latent dimension ", dim2),
     main = "eSVD embedding and density\n(Curved Gaussian)",
     xlim = range(x), ylim = range(y))

for(ll in plotting_order){
  target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx == ll)]
  idx <- which(cluster_labels %in% target_indices)
  points(x = x[idx], y = y[idx], pch = 16,
         col = col_vec2_esvd[target_indices[1]])
}

quantile_vec <- rep(0.95, length(plotting_order_esvd))
quantile_vec[4] <- 0.99
quantile_vec[2] <- 0.94

for(ll in plotting_order){
  target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx == ll)]
  idx <- which(cluster_labels %in% target_indices)
  kde_est <- MASS::kde2d(x[idx], y[idx], n = 500, lims = c(range(x), range(y)))
  col_val <- col_info_esvd[target_indices[1], "col_code"]
  quant_level <- quantile(kde_est$z, probs = quantile_vec[ll])
  contour(kde_est, add = T, drawlabels = F, lwd = 7, levels = quant_level, col = "white")
  contour(kde_est, add = T, drawlabels = F, lwd = 4, levels = quant_level, col = col_val)
  contour(kde_est, add = T, drawlabels = F, lwd = 2, levels = quant_level,  lty = 3, col = "black")
}
graphics.off()
