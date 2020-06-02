var <- ls()
rm(list = var[var != "suffix"])
load(paste0("../results/step7_figures", suffix, ".RData"))

esvd_curves_prepared <- lapply(esvd_curves_short$curves, function(curve){
  prepare_trajectory(curve, target_length = 11)
})

esvd_tube_list <- lapply(1:length(esvd_curves_prepared), function(x){
  s_mat <- esvd_curves_prepared[[x]][,1:3]
  eSVD::construct_3d_tube(s_mat, radius = esvd_width)
})

svd_tube_list <- lapply(1:length(svd_curves_short$curves), function(x){
  s_mat <- svd_curves_short$curves[[x]]$s[svd_curves_short$curves[[x]]$ord,][,1:3]
  eSVD::construct_3d_tube(s_mat, radius = svd_width)
})

###########################################

#angle_matrix
angle_matrix <- matrix(c(90,225, 180,225, 270,225), byrow = T, ncol = 2)
spacing <- 4.5; tmp <- spacing*c(-.5,.5)
bound_matrix <- list(list(xlim = -2+tmp, ylim = 0+tmp, zlim = 1+tmp),
                     list(xlim = -2.5+tmp, ylim = 0.75+tmp, zlim = 1+tmp),
                     list(xlim = -3+tmp, ylim = 0.75+tmp, zlim = 1+tmp))
tmp <- color_func(0.07); tmp[6] <- rgb(100/255, 100/255, 100/255, 0.05)
col_vec3_esvd <- tmp[num_order_vec_esvd]

for(kk in 1:nrow(angle_matrix)){
  print(kk)
  png(paste0("../../esvd_results/figure/main/eSVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0,0,4,0))
  eSVD:::slingshot_3dplot(esvd_embedding$u_mat[,1:3], cluster_labels,
                          bg_col_vec = col_vec3_esvd, bg_cex = 0.8,
                          breaks = seq(0.5, 13.5, by = 1),
                          cluster_center = cluster_center_esvd,
                          center_col_vec = col_vec_esvd,
                          center_labels = 1:13,
                          curves = NA,
                          pch = 16, main = "eSVD embedding and trajectories:\nCurved Gaussian",
                          xlab = "Latent dimension 1", ylab = "Latent dimension 2",
                          zlab = "Latent dimension 3",
                          theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                          xlim = bound_matrix[[kk]]$xlim,
                          ylim = bound_matrix[[kk]]$ylim,
                          zlim = bound_matrix[[kk]]$zlim)

  curves <- esvd_curves_prepared
  col_vec_short <- color_func(1)[c(1,4)]

  for(i in 1:length(curves)){
    plot3D::lines3D(x = curves[[i]][, 1],
                    y = curves[[i]][, 2],
                    z = curves[[i]][, 3],
                    add = T, colkey = F, col = "black", lwd = 12)
    plot3D::lines3D(x = curves[[i]][, 1],
                    y = curves[[i]][, 2],
                    z = curves[[i]][, 3],
                    add = T, colkey = F, col = col_vec_short[i], lwd = 12)
  }
  graphics.off()
}

# 3D plots with tubes

for(kk in 1:nrow(angle_matrix)){
  print(kk)
  png(paste0("../../esvd_results/figure/main/eSVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], "_tube.png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0,0,4,0))
  eSVD:::slingshot_3dplot(esvd_embedding$u_mat[,1:3], cluster_labels,
                          bg_col_vec = col_vec3_esvd, bg_cex = 0.8,
                          breaks = seq(0.5, 13.5, by = 1),
                          cluster_center = cluster_center_esvd,
                          center_col_vec = col_vec_esvd,
                          center_labels = 1:13,
                          curves = NA,
                          pch = 16, main = "eSVD embedding and trajectories:\nCurved Gaussian",
                          xlab = "Latent dimension 1", ylab = "Latent dimension 2",
                          zlab = "Latent dimension 3",
                          theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                          xlim = bound_matrix[[kk]]$xlim,
                          ylim = bound_matrix[[kk]]$ylim,
                          zlim = bound_matrix[[kk]]$zlim)

  curves <- esvd_curves_prepared
  col_vec_short <- color_func(0.9)[c(1,4)]

  for(i in 1:length(curves)){
    plot3D::lines3D(x = curves[[i]][, 1],
                    y = curves[[i]][, 2],
                    z = curves[[i]][, 3],
                    add = T, colkey = F, col = "black", lwd = 6)
    plot3D::lines3D(x = curves[[i]][, 1],
                    y = curves[[i]][, 2],
                    z = curves[[i]][, 3],
                    add = T, colkey = F, col = col_vec_short[i], lwd = 6)

    col_mat <- esvd_tube_list[[i]]$z_mat

    plot3D::surf3D(esvd_tube_list[[i]]$x_mat,
                   esvd_tube_list[[i]]$y_mat,
                   esvd_tube_list[[i]]$z_mat, add = T,
                   colvar = col_mat,
                   col = colorRampPalette(c("white", col_vec_short[i]))(150)[51:150],
                   breaks = seq(min(col_mat), max(col_mat), length.out = 101),
                   colkey = F)
  }
  graphics.off()
}


###############################

#angle_matrix
angle_matrix <- matrix(c(0,135, 90,135, 225,45), byrow = T, ncol = 2)
spacing <- 2.25; tmp <- spacing*c(-.5,.5)
bound_matrix <- list(list(xlim = -1.5+tmp, ylim = -0.5+tmp, zlim = 0+tmp),
                     list(xlim = -1.5+tmp, ylim = -0.5+tmp, zlim = 0+tmp),
                     list(xlim = -1.5+tmp, ylim = -0.5+tmp, zlim = 0+tmp))
col_vec3_svd <- color_func(0.09)[num_order_vec_svd]

for(kk in 1:nrow(angle_matrix)){
  print(kk)
  png(paste0("../../esvd_results/figure/main/SVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0,0,4,0))
  eSVD:::slingshot_3dplot(svd_embedding[,1:3], cluster_labels,
                          bg_col_vec = col_vec3_svd, bg_cex = 0.8,
                          breaks = seq(0.5, 13.5, by = 1),
                          cluster_center = cluster_center_svd,
                          center_col_vec = col_vec_svd,
                          center_labels = 1:13,
                          curves = NA,
                          pch = 16, main = "SVD embedding and trajectories:\nConstant-variance Gaussian",
                          xlab = "Latent dimension 1", ylab = "Latent dimension 2",
                          zlab = "Latent dimension 3",
                          theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                          xlim = bound_matrix[[kk]]$xlim,
                          ylim = bound_matrix[[kk]]$ylim,
                          zlim = bound_matrix[[kk]]$zlim)

  curves <- svd_curves_short$curves

  for(i in 1:length(curves)){
    ord <- curves[[i]]$ord
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = "black", lwd = 6)
  }
  graphics.off()
}

# 3D plots with tubes
col_vec_short <- color_func(0.5)[c(4)]

for(kk in 1:nrow(angle_matrix)){
  print(kk)
  png(paste0("../../esvd_results/figure/main/SVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], "_tube.png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0,0,4,0))
  eSVD:::slingshot_3dplot(svd_embedding[,1:3], cluster_labels,
                          bg_col_vec = col_vec3_svd, bg_cex = 0.8,
                          breaks = seq(0.5, 13.5, by = 1),
                          cluster_center = cluster_center_svd,
                          center_col_vec = col_vec_svd,
                          center_labels = 1:13,
                          curves = NA,
                          pch = 16, main = "SVD embedding and trajectories:\nConstant-variance Gaussian",
                          xlab = "Latent dimension 1", ylab = "Latent dimension 2",
                          zlab = "Latent dimension 3",
                          theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                          xlim = bound_matrix[[kk]]$xlim,
                          ylim = bound_matrix[[kk]]$ylim,
                          zlim = bound_matrix[[kk]]$zlim)

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
                 col = (colorRampPalette(c(col_vec_short[1], "white"))(110))[1:100],
                 breaks = seq(min(col_mat), max(col_mat), length.out = 101),
                 colkey = F)
  graphics.off()
}

