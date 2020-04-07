rm(list=ls())
load("../results/step5_trajectory.RData")
esvd_curves$lineages
esvd_sd_val$sd_val
svd_sd_val$sd_val

# start with 2D plots

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

#####################


num_order_vec_esvd <- c(5, rep(3,2), c(6,1,1,1,6,4), rep(2,2),  rep(5,2))
col_vec_esvd <- color_func(1)[num_order_vec_esvd]
col_vec2_esvd <- color_func(0.5)[num_order_vec_esvd]
col_name_esvd <- c("orange", rep("bluish green", 2), c("gray", "yellow", "yellow", "yellow", "gray", "blue"), rep("skyblue", 2), rep("orange", 2))
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
plotting_order <- c(3,4,5,2,6,1)

cluster_center_esvd <- .compute_cluster_center(esvd_embedding$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))

######

dim1 <- 1; dim2 <- 3

png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup6_main_esvd_2dplots_1and3.png"),
    height = 1500, width = 1500, res = 300,
    units = "px")
plot(NA, xlim = range(esvd_embedding$u_mat[,dim1]), ylim = range(esvd_embedding$u_mat[,dim2]),
     asp = T, xlab = paste0("Latent dimension ", dim1), ylab = paste0("Latent dimension ", dim2),
     main = "eSVD embedding and trajectories\n(Curved Gaussian)")

for(ll in plotting_order) {
  target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx == ll)]
  idx <- which(cluster_labels %in% target_indices)
  points(x = esvd_embedding$u_mat[idx,dim1], y = esvd_embedding$u_mat[idx,dim2], pch = 16,
         col = col_vec2_esvd[target_indices[1]])
}

for(ll in 1:nrow(cluster_center_esvd)){
  points(cluster_center_esvd[ll,dim1], cluster_center_esvd[ll,dim2], pch = 16, cex = 2, col = "black")
  points(cluster_center_esvd[ll,dim1], cluster_center_esvd[ll,dim2], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
}


curves <- esvd_curves$curves
for(ll in 1:length(curves)) {
  ord <- curves[[ll]]$ord
  lines(x = curves[[ll]]$s[ord, dim1], y = curves[[ll]]$s[ord, dim2], col = "white", lwd = 8)
  lines(x = curves[[ll]]$s[ord, dim1], y = curves[[ll]]$s[ord, dim2], col = col_vec_short[ll], lwd = 5)
  lines(x = curves[[ll]]$s[ord, dim1], y = curves[[ll]]$s[ord, dim2], col = "black",
        lty = 3, lwd = 2)
}

graphics.off()

######

png(paste0("../../esvd_results/figure/experiment/Revision_writeup6_main_esvd_density.png"),
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

quantile_vec <- rep(0.95, length(plotting_order))
quantile_vec[4] <- 0.99
quantile_vec[2] <- 0.94

for(ll in plotting_order){
  target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx == ll)]
  idx <- which(cluster_labels %in% target_indices)
  kde_est <- MASS::kde2d(x[idx], y[idx], n = 500, lims = c(range(x), range(y)))
  col_val <- col_info_esvd[target_indices[1], "col_code"]
  quant_level <- quantile(kde_est$z, probs = quantile_vec[ll])
  contour(kde_est, add = T, drawlabels = F, lwd = 7, levels = quant_level, col = "white")
  #contour(kde_est, add = T, drawlabels = F, lwd = 2.5, levels = quant_level, col = "black")
  contour(kde_est, add = T, drawlabels = F, lwd = 4, levels = quant_level, col = col_val)
  contour(kde_est, add = T, drawlabels = F, lwd = 2, levels = quant_level,  lty = 3, col = "black")
}
graphics.off()
