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

############

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

combn_mat <- combn(3,2)
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup6_main_esvd_2dplots_", k, ".png"),
      height = 1750, width = 1750, res = 300,
      units = "px")
  plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding and trajectories\n(Curved Gaussian)")

  for(ll in plotting_order) {
    target_indices <- col_info_esvd$idx[which(col_info_esvd$factor_idx == ll)]
    idx <- which(cluster_labels %in% target_indices)
    points(x = esvd_embedding$u_mat[idx,i], y = esvd_embedding$u_mat[idx,j], pch = 16,
           col = col_info_esvd$col_code[target_indices[1]])
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

  graphics.off()
}

#############################

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
plotting_order <- c(2,3,1,4)

cluster_center_svd <- .compute_cluster_center(svd_embedding[,1:3], .construct_cluster_matrix(cluster_labels))

combn_mat <- combn(3,2)
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  png(filename = paste0("../../esvd_results/figure/experiment/Revision_writeup6_main_svd_2dplots_", k, ".png"),
      height = 1750, width = 1750, res = 300,
      units = "px")
  plot(NA, xlim = range(svd_embedding[,i]), ylim = range(svd_embedding[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "SVD embedding and trajectories\n(Constant-variance Gaussian)")

  for(ll in plotting_order) {
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

#####################################################################


#angle_matrix
angle_matrix = matrix(c(45,225, 270,135, 225,360), byrow = T, ncol = 2)
bound_matrix = list(list(xlim = c(0, 3), ylim = c(-2, 2), zlim = c(-1.5, 1)),
                    list(xlim = c(0, 3), ylim = c(-2, 2), zlim = c(-1.5, 1)),
                    list(xlim = c(0, 3), ylim = c(-2, 2), zlim = c(-1.5, 1)))
col_vec3_esvd <- color_func(0.09)[num_order_vec_esvd]

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
angle_matrix = matrix(c(45,225, 270,135, 225,360), byrow = T, ncol = 2)
bound_matrix = list(list(xlim = c(-8, 0), ylim = c(-6, 4), zlim = c(-1, 20)),
                    list(xlim = c(-8, 0), ylim = c(-6, 4), zlim = c(-1, 20)),
                    list(xlim = c(-8, 0), ylim = c(-6, 4), zlim = c(-1, 20)))
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
col_vec_short <- color_func(0.9)[c(4)]

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
                 col = colorRampPalette(c("white", col_vec_short[1]))(100),
                 breaks = seq(min(col_mat), max(col_mat), length.out = 101),
                 colkey = F)
  graphics.off()
}

