rm(list=ls())
load("../results/step5_trajectory_original.RData")

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(100/255, 100/255, 100/255, alpha)) #gray
}
color_name_vec <- c("yellow", "skyblue", "bluish green", "blue", "orange", "gray")


cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

num_order_vec_esvd <- c(5, rep(3,2), c(6,1,1,4,4,4), rep(2,2),  rep(5,2))
col_vec_esvd <- color_func(1)[num_order_vec_esvd]
col_vec2_esvd <- color_func(0.5)[num_order_vec_esvd]
col_name_esvd <- color_name_vec[num_order_vec_esvd]
order_vec_esvd <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
col_info_esvd <- data.frame(name = levels(cell_type_vec),
                            idx = sort(unique(cluster_labels)),
                            order = order_vec_esvd,
                            col_name = col_name_esvd,
                            col_code = col_vec_esvd)
col_info_esvd$factor_idx <- as.numeric(as.factor(col_info_esvd$col_name))
col_info_esvd[,c(5,6)] <- col_info_esvd[,c(6,5)]
colnames(col_info_esvd)[c(5,6)] <- colnames(col_info_esvd)[c(6,5)]
col_vec_short <- color_func(0.9)[c(1,4)]
cluster_center_esvd <- eSVD:::.compute_cluster_center(esvd_embedding$u_mat[,1:3], eSVD:::.construct_cluster_matrix(cluster_labels))


num_order_vec_svd <- c(5, rep(3,2), c(1,1,1,1,1,1), rep(2,2),  rep(5,2))
col_vec_svd <- color_func(1)[num_order_vec_svd]
col_vec2_svd <- color_func(0.5)[num_order_vec_svd]
col_name_svd <- color_name_vec[num_order_vec_svd]
order_vec_svd <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
col_info_svd <- data.frame(name = levels(cell_type_vec),
                           idx = sort(unique(cluster_labels)),
                           order = order_vec_svd,
                           col_name = col_name_svd,
                           col_code = col_vec_svd)
col_info_svd$factor_idx <- as.numeric(as.factor(col_info_svd$col_name))
col_info_svd[,c(5,6)] <- col_info_svd[,c(6,5)]
colnames(col_info_svd)[c(5,6)] <- colnames(col_info_svd)[c(6,5)]
cluster_center_svd <- eSVD:::.compute_cluster_center(svd_embedding[,1:3], eSVD:::.construct_cluster_matrix(cluster_labels))

bound_matrix_esvd <- list(list(xlim = c(-5, 1), ylim = c(-2.5, 3.5), zlim = c(-2, 4)),
                          list(xlim = c(-5, 1), ylim = c(-2.5, 3.5), zlim = c(-2, 4)),
                          list(xlim = c(-5, 1), ylim = c(-2.5, 3.5), zlim = c(-2, 4)))
col_vec3_esvd <- color_func(0.09)[num_order_vec_esvd]

bound_matrix_svd <- list(list(xlim = c(-2.5, 0.5), ylim = c(-2, 1), zlim = c(-1.5, 1.5)),
                         list(xlim = c(-2.5, 0.5), ylim = c(-2, 1), zlim = c(-1.5, 1.5)),
                         list(xlim = c(-2.5, 0.5), ylim = c(-2, 1), zlim = c(-1.5, 1.5)))
col_vec3_svd <- color_func(0.09)[num_order_vec_svd]

esvd_width <- max(sapply(1:2, function(i){
  quantile(apply(esvd_sd_res$mat_list[[i]], 1, function(x){quantile(x, probs = 0.95)}), probs = 0.95)
}))
svd_width <- max(sapply(1:2, function(i){
  quantile(apply(svd_sd_val$mat_list[[i]], 1, function(x){quantile(x, probs = 0.95)})[200:900], probs = 0.95)
}))

esvd_tube_list <- lapply(1:length(esvd_curves_short$curves), function(x){
  s_mat <- esvd_curves_short$curves[[x]]$s[esvd_curves_short$curves[[x]]$ord,][,1:3]
  eSVD::construct_3d_tube(s_mat, radius = esvd_width)
})
svd_tube_list <- lapply(1:length(svd_curves_short$curves), function(x){
  s_mat <- svd_curves_short$curves[[x]]$s[svd_curves_short$curves[[x]]$ord,][,1:3]
  eSVD::construct_3d_tube(s_mat, radius = svd_width)
})

########################

angle_matrix <- as.matrix(expand.grid(seq(0, 360, length.out = 9), seq(0, 360, length.out = 9)))

for(kk in 1:nrow(angle_matrix)){
  png(paste0("../../esvd_results/figure/3dplots/eSVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0,0,4,0))
  eSVD:::slingshot_3dplot(esvd_embedding$u_mat[,1:3], cluster_labels,
                          bg_col_vec = col_vec3_esvd, bg_cex = 0.8,
                          breaks = seq(0.5, 13.5, by = 1),
                          cluster_center = cluster_center_esvd,
                          center_col_vec = col_vec_esvd,
                          center_labels = 1:13,
                          curves = NA,
                          pch = 16, main = "eSVD embedding and trajectories\n(Curved Gaussian)",
                          xlab = "Latent dimension 1", ylab = "Latent dimension 2",
                          zlab = "Latent dimension 3",
                          theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                          xlim = bound_matrix_esvd[[1]]$xlim,
                          ylim = bound_matrix_esvd[[1]]$ylim,
                          zlim = bound_matrix_esvd[[1]]$zlim)

  curves <- esvd_curves_short$curves
  col_vec_short <- color_func(0.9)[c(1,4)]

  for(i in 1:length(curves)){
    ord <- curves[[i]]$ord
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = "black", lwd = 6)
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = col_vec_short[i], lwd = 6)
  }
  graphics.off()

  #############

  png(paste0("../../esvd_results/figure/3dplots/eSVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], "_tube.png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0,0,4,0))
  eSVD:::slingshot_3dplot(esvd_embedding$u_mat[,1:3], cluster_labels,
                          bg_col_vec = col_vec3_esvd, bg_cex = 0.8,
                          breaks = seq(0.5, 13.5, by = 1),
                          cluster_center = cluster_center_esvd,
                          center_col_vec = col_vec_esvd,
                          center_labels = 1:13,
                          curves = NA,
                          pch = 16, main = "eSVD embedding and trajectories\n(Curved Gaussian)",
                          xlab = "Latent dimension 1", ylab = "Latent dimension 2",
                          zlab = "Latent dimension 3",
                          theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                          xlim = bound_matrix_esvd[[1]]$xlim,
                          ylim = bound_matrix_esvd[[1]]$ylim,
                          zlim = bound_matrix_esvd[[1]]$zlim)

  curves <- esvd_curves_short$curves
  col_vec_short <- color_func(0.9)[c(1,4)]

  for(i in 1:length(curves)){
    ord <- curves[[i]]$ord
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = "black", lwd = 6)
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = col_vec_short[i], lwd = 6)

    col_mat <- esvd_tube_list[[i]]$z_mat

    plot3D::surf3D(esvd_tube_list[[i]]$x_mat,
                   esvd_tube_list[[i]]$y_mat,
                   esvd_tube_list[[i]]$z_mat, add = T,
                   colvar = col_mat,
                   col = colorRampPalette(c("white", col_vec_short[i]))(100),
                   breaks = seq(min(col_mat), max(col_mat), length.out = 101),
                   colkey = F)
  }
  graphics.off()

  ##############

  png(paste0("../../esvd_results/figure/3dplots/SVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0,0,4,0))
  eSVD:::slingshot_3dplot(svd_embedding[,1:3], cluster_labels,
                          bg_col_vec = col_vec3_svd, bg_cex = 0.8,
                          breaks = seq(0.5, 13.5, by = 1),
                          cluster_center = cluster_center_svd,
                          center_col_vec = col_vec_svd,
                          center_labels = 1:13,
                          curves = NA,
                          pch = 16, main = "SVD embedding and trajectories\n(Constant-variance Gaussian)",
                          xlab = "Latent dimension 1", ylab = "Latent dimension 2",
                          zlab = "Latent dimension 3",
                          theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                          xlim = bound_matrix_svd[[1]]$xlim,
                          ylim = bound_matrix_svd[[1]]$ylim,
                          zlim = bound_matrix_svd[[1]]$zlim)

  curves <- svd_curves_short$curves

  for(i in 1:length(curves)){
    ord <- curves[[i]]$ord
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = "black", lwd = 6)
  }
  graphics.off()

  #####

  png(paste0("../../esvd_results/figure/3dplots/SVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], "_tube.png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0,0,4,0))
  eSVD:::slingshot_3dplot(svd_embedding[,1:3], cluster_labels,
                          bg_col_vec = col_vec3_svd, bg_cex = 0.8,
                          breaks = seq(0.5, 13.5, by = 1),
                          cluster_center = cluster_center_svd,
                          center_col_vec = col_vec_svd,
                          center_labels = 1:13,
                          curves = NA,
                          pch = 16, main = "SVD embedding and trajectories\n(Constant-variance Gaussian)",
                          xlab = "Latent dimension 1", ylab = "Latent dimension 2",
                          zlab = "Latent dimension 3",
                          theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                          xlim = bound_matrix_svd[[1]]$xlim,
                          ylim = bound_matrix_svd[[1]]$ylim,
                          zlim = bound_matrix_svd[[1]]$zlim)

  curves <- svd_curves_short$curves

  for(i in 1:length(curves)){
    ord <- curves[[i]]$ord
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = "black", lwd = 6)
  }

  i <- 2
  col_mat <- svd_tube_list[[i]]$z_mat

  plot3D::surf3D(svd_tube_list[[i]]$x_mat,
                 svd_tube_list[[i]]$y_mat,
                 svd_tube_list[[i]]$z_mat, add = T,
                 colvar = col_mat,
                 col = colorRampPalette(c("white", col_vec_short[2]))(100),
                 breaks = seq(min(col_mat), max(col_mat), length.out = 101),
                 colkey = F)
  graphics.off()

}

