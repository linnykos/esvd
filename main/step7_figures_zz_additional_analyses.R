var <- ls()
rm(list = var[var != "suffix"])
load(paste0("../results/step7_figures", suffix, ".RData"))
cluster_center_zinbwave <- eSVD::compute_cluster_center(zinbwave_embedding[,1:3], .construct_cluster_matrix(cluster_labels))

## zinbwave plots
grDevices::png(filename = paste0("../../esvd_results/figure/main/zinbwave_2dplots.png"),
               height = 830, width = 2300, res = 300,
               units = "px")
graphics::par(mfrow = c(1,3), mar = c(4,4,4,1))

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  graphics::plot(NA, xlim = range(zinbwave_embedding[,i]), ylim = range(zinbwave_embedding[,j]),
                 asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
                 main = ifelse(k == 2, "ZINB-WaVE embedding", ""))

  for(ll in plotting_order_esvd) {
    target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx %in% ll)]
    idx <- which(cluster_labels %in% target_indices)
    graphics::points(x = zinbwave_embedding[idx,i], y = zinbwave_embedding[idx,j], pch = 16,
                     col = col_vec2_esvd[cluster_labels[idx]])
  }

  for(ll in 1:nrow(cluster_center_zinbwave)){
    graphics::points(cluster_center_zinbwave[ll,i], cluster_center_zinbwave[ll,j], pch = 16, cex = 2, col = "black")
    graphics::points(cluster_center_zinbwave[ll,i], cluster_center_zinbwave[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
  }
}
grDevices::graphics.off()

################################

## UMAP plots
cluster_center_umap <- eSVD::compute_cluster_center(umap_all$layout, .construct_cluster_matrix(cluster_labels))

grDevices::png(filename = paste0("../../esvd_results/figure/main/umap.png"),
               height = 1500, width = 2500, res = 300,
               units = "px")
graphics::par(mfrow = c(1,2))
graphics::plot(NA, asp = T, pch = 16, xlim = range(umap_all$layout[,1]), ylim = range(umap_all$layout[,2]),
               xlab = "UMAP dimension 1", ylab = "UMAP dimension 2", main = "UMAP embedding and centers")
for(ll in plotting_order_esvd){
  target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx %in% ll)]
  idx <- which(cluster_labels %in% target_indices)
  graphics::points(x = umap_all$layout[idx,1], y = umap_all$layout[idx,2], pch = 16,
                   col = col_vec2_esvd[cluster_labels[idx]])
}

for(ll in 1:nrow(cluster_center_umap)){
  graphics::points(cluster_center_umap[ll,1], cluster_center_umap[ll,2], pch = 16, cex = 2.25, col = "black")
  graphics::points(cluster_center_umap[ll,1], cluster_center_umap[ll,2], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
}

graphics::plot(NA, asp = T, pch = 16, xlim = range(umap_all$layout[,1]), ylim = range(umap_all$layout[,2]),
               xlab = "UMAP dimension 1", ylab = "UMAP dimension 2", main = "UMAP embedding and densities")
for(ll in plotting_order_esvd){
  target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx %in% ll)]
  idx <- which(cluster_labels %in% target_indices)
  graphics::points(x = umap_all$layout[idx,1], y = umap_all$layout[idx,2], pch = 16,
                   col = col_vec2_esvd[cluster_labels[idx]])
}

quantile_vec <- rep(0.95, 6)
x <- umap_all$layout[,1]; y <- umap_all$layout[,2]

for(ll in 1:6){
  target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx %in% ll)]
  idx <- which(cluster_labels %in% target_indices)
  kde_est <- MASS::kde2d(x[idx], y[idx], n = 500, lims = c(range(x), range(y)))
  col_val <- col_info_esvd[target_indices[1], "col_code"]
  quant_level <- stats::quantile(kde_est$z, probs = quantile_vec[ll])
  graphics::contour(kde_est, add = T, drawlabels = F, lwd = 7, levels = quant_level, col = "white")
  graphics::contour(kde_est, add = T, drawlabels = F, lwd = 4, levels = quant_level, col = col_val)
  graphics::contour(kde_est, add = T, drawlabels = F, lwd = 2, levels = quant_level,  lty = 3, col = "black")
}

grDevices::graphics.off()

########################################

nat_mat_list <- lapply(1:length(esvd_missing_list2[[esvd_angle_res2$idx]]), function(i){
  esvd_missing_list2[[esvd_angle_res2$idx]][[i]]$u_mat %*% t(esvd_missing_list2[[esvd_angle_res2$idx]][[i]]$v_mat)
})

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

png(filename = paste0("../../esvd_results/figure/main/esvd_training_testing_nb.png"),
    height = 1500, width = 2500, res = 300,
    units = "px")
par(mfrow = c(1,2))
eSVD::plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                       missing_idx_list = training_idx_list,
                                       family = "neg_binom", cex.lab = 1.25,
                                       scalar = paramMat_esvd2[esvd_angle_res2$idx, "scalar"],
                                       main = "eSVD embedding's diagnostic:\nNegative binomial,\ntraining set",
                                       max_points = 1e6)


eSVD::plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                       missing_idx_list = missing_idx_list,
                                       family = "neg_binom", cex.lab = 1.25,
                                       scalar = paramMat_esvd2[esvd_angle_res2$idx, "scalar"],
                                       main = "eSVD embedding's diagnostic:\nNegative binomial,\ntesting set")
graphics.off()
