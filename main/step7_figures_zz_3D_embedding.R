var <- ls()
rm(list = var[var != "suffix"])
load(paste0("../results/step7_figures", suffix, ".RData"))


###########################################

#angle_matrix
angle_matrix <- matrix(c(0,45, 180,0, 270,45), byrow = T, ncol = 2)
spacing <- 3
bound_matrix <- list(list(xlim = 0.5+c(0,spacing), ylim = -2+c(0,spacing), zlim = -1.5+c(0,spacing)),
                     list(xlim = 0.5+c(0,spacing), ylim = -2+c(0,spacing), zlim = -2+c(0,spacing)),
                     list(xlim = 0.5+c(0,spacing), ylim = -1+c(0,spacing), zlim = -1.5+c(0,spacing)))
tmp <- color_func(0.07); tmp[6] <- rgb(100/255, 100/255, 100/255, 0.05)
col_vec3_esvd <- tmp[num_order_vec_esvd]

for(kk in 1:nrow(angle_matrix)){
  print(kk)
  png(paste0("../../esvd_results/figure/experiment/eSVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], ".png"),
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
                          xlim = bound_matrix[[kk]]$xlim,
                          ylim = bound_matrix[[kk]]$ylim,
                          zlim = bound_matrix[[kk]]$zlim)

  curves <- esvd_curves$curves
  col_vec_short <- color_func(1)[c(1,4)]

  for(i in 1:length(curves)){
    ord <- curves[[i]]$ord
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = "black", lwd = 12)
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = col_vec_short[i], lwd = 12)
  }
  graphics.off()
}

# 3D plots with tubes

esvd_tube_list <- lapply(1:length(esvd_curves$curves), function(x){
  s_mat <- esvd_curves$curves[[x]]$s[esvd_curves$curves[[x]]$ord,]
  eSVD::construct_3d_tube(s_mat, radius = esvd_sd_val$sd_val/2)
})

for(kk in 1:nrow(angle_matrix)){
  print(kk)
  png(paste0("../../esvd_results/figure/experiment/eSVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], "_tube.png"),
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
                          xlim = bound_matrix[[kk]]$xlim,
                          ylim = bound_matrix[[kk]]$ylim,
                          zlim = bound_matrix[[kk]]$zlim)

  curves <- esvd_curves$curves
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
}


###############################

#angle_matrix
angle_matrix <- matrix(c(45,0, 135,0,  225,0), byrow = T, ncol = 2)
spacing <- 8
bound_matrix <- list(list(xlim = -8+c(0,spacing), ylim = -6+c(0,spacing), zlim = 0+c(0,spacing)),
                     list(xlim = -7+c(0,spacing), ylim = -5+c(0,spacing), zlim = -1+c(0,spacing)),
                     list(xlim = -8+c(0,spacing), ylim = -3+c(0,spacing), zlim = 0+c(0,spacing)))
col_vec3_svd <- color_func(0.09)[num_order_vec_svd]

for(kk in 1:nrow(angle_matrix)){
  print(kk)
  png(paste0("../../esvd_results/figure/experiment/SVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], ".png"),
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
                          xlim = bound_matrix[[kk]]$xlim,
                          ylim = bound_matrix[[kk]]$ylim,
                          zlim = bound_matrix[[kk]]$zlim)

  curves <- svd_curves$curves

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

svd_tube_list <- lapply(1:length(svd_curves$curves), function(x){
  s_mat <- svd_curves$curves[[x]]$s[svd_curves$curves[[x]]$ord,]
  eSVD::construct_3d_tube(s_mat, radius = svd_sd_val$sd_val/2)
})
col_vec_short <- color_func(0.5)[c(4)]

for(kk in 1:nrow(angle_matrix)){
  print(kk)
  png(paste0("../../esvd_results/figure/experiment/SVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], "_tube.png"),
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
                          xlim = bound_matrix[[kk]]$xlim,
                          ylim = bound_matrix[[kk]]$ylim,
                          zlim = bound_matrix[[kk]]$zlim)

  curves <- svd_curves$curves

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
                 breaks = c(seq(min(col_mat), 8, length.out = 100), max(col_mat)),
                 colkey = F)
  graphics.off()
}

