var <- ls()
rm(list = var[var != "suffix"])
load(paste0("../results/step6_figures", suffix, ".RData"))

dim1 <- 2; dim2 <- 3

png(paste0("../../esvd_results/figure/main/svd_density.png"),
    height = 1300, width = 1500, res = 300, units = "px")
x <- svd_embedding[,dim1]
y <- svd_embedding[,dim2]
graphics::plot(NA,  asp = T, xlab = paste0("Latent dimension ", dim1), ylab = paste0("Latent dimension ", dim2),
     main = "SVD embedding\n(Constant-variance Gaussian)",
     xlim = range(x), ylim = range(y), axes = F)
graphics::points(x = x, y = y, col = col_vec2_svd[cluster_labels], pch = 16)
graphics::axis(1); graphics::axis(2)

for(i in 1:4){
  clust_idx <- which(col_info_svd[,"factor_idx"] %in% i)
  idx <- which(cluster_labels %in% clust_idx)
  kde_est <- MASS::kde2d(svd_embedding[idx,dim1], svd_embedding[idx,dim2], n = 500,
                         lims = c(range(x), range(y)))
  col_val <- col_info_svd[clust_idx[1], "col_code"]
  quant_level <- stats::quantile(kde_est$z, probs = 0.925)
  graphics::contour(kde_est, add = T, drawlabels = F, lwd = 7, levels = quant_level, col = "white")
  graphics::contour(kde_est, add = T, drawlabels = F, lwd = 4, levels = quant_level, col = col_val)
  graphics::contour(kde_est, add = T, drawlabels = F, lwd = 2, levels = quant_level,  lty = 3, col = "black")
}
graphics.off()


###########################

# esvd density

dim1 <- 2
dim2 <- 3
png(paste0("../../esvd_results/figure/main/esvd_density.png"),
    height = 1500, width = 1500, res = 300, units = "px")
x <- esvd_embedding$u_mat[,dim1]
y <- esvd_embedding$u_mat[,dim2]
graphics::plot(NA,  asp = T, xlab = paste0("Latent dimension ", dim1), ylab = paste0("Latent dimension ", dim2),
     main = "eSVD embedding and density\n(Curved Gaussian)",
     xlim = range(x), ylim = range(y))

for(ll in plotting_order_esvd){
  target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx == ll)]
  idx <- which(cluster_labels %in% target_indices)
  graphics::points(x = x[idx], y = y[idx], pch = 16,
         col = col_vec2_esvd[target_indices[1]])
}

quantile_vec <- rep(0.99, length(plotting_order_esvd))
quantile_vec[3] <- 0.975

for(ll in c(3,6,1)){
  target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx == ll)]
  idx <- which(cluster_labels %in% target_indices)
  kde_est <- MASS::kde2d(x[idx], y[idx], n = 500, lims = c(range(x), range(y)))
  col_val <- col_info_esvd[target_indices[1], "col_code"]
  quant_level <- stats::quantile(kde_est$z, probs = quantile_vec[ll])
  graphics::contour(kde_est, add = T, drawlabels = F, lwd = 7, levels = quant_level, col = "white")
  graphics::contour(kde_est, add = T, drawlabels = F, lwd = 4, levels = quant_level, col = col_val)
  graphics::contour(kde_est, add = T, drawlabels = F, lwd = 2, levels = quant_level,  lty = 3, col = "black")
}
graphics.off()
