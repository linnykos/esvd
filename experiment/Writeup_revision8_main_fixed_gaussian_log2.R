rm(list=ls())
load("../results/step3_scalar_heuristic_cg_hvg.RData")

# first fit the trajectories
cell_type_vec <- as.character(marques$cell.info$cell.type)
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

upscale_factor <- 0.5
reduction_percentage <- 0.3
p <- 3

set.seed(10)
svd_curves <- slingshot(svd_embedding[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                        cluster_group_list = cluster_group_list,
                        verbose = T, upscale_factor = upscale_factor, reduction_percentage = reduction_percentage, squared = F)


###############

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

cluster_center_svd <- .compute_cluster_center(svd_embedding[,1:3], .construct_cluster_matrix(cluster_labels))
combn_mat <- combn(3,2)

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision8_svd_2dplots_", k, ".png"),
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
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 5)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black", lwd = 3)
  }

  graphics.off()
}



