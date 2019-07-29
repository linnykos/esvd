# load(paste0("../results/step5_clustering", suffix, ".RData"))
load(paste0("../results/old_results/step5_clustering_spca.RData"))

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(210/255, 198/255, 36/255, alpha)) #darker yellow
}
cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)

######################################

# first do our embedding

cluster_center <- .compute_cluster_center(res_our$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))
custom_cluster_group_list <- list(13, 12, 1, c(10,11), c(2,3), c(4:7), c(8,9))

num_order_vec <- c(5, rep(3,2), 3, rep(1,3), rep(4,2), rep(2,2),  rep(5,2))
col_vec <- color_func(1)[num_order_vec]
col_vec2 <- color_func(0.1)[num_order_vec]
col_name <- c("orange", rep("bluish green", 3), rep("yellow", 3), rep("blue", 2), rep("skyblue", 2), rep("orange", 2))
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                       level = sapply(1:13, function(y){which(sapply(custom_cluster_group_list, function(x){y %in% x}))}),
                       col_name = col_name,
                       col_code = col_vec)
col_info

#angle_matrix
angle_matrix = matrix(c(45,225, 270,135, 225,360), byrow = T, ncol = 2)
bound_matrix = list(list(xlim = c(-2.45, 0), ylim = c(-1.5,1), zlim = c(-1,0.75)),
                    list(xlim = c(-2.45, 0), ylim = c(-1.5,1), zlim = c(-1,0.75)),
                    list(xlim = c(-2.45, 0), ylim = c(-1.5,1), zlim = c(-1,0.75)))

for(kk in 1:nrow(angle_matrix)){
  print(kk)
  png(paste0("../figure/main/eSVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0,0,4,0))
  slingshot_3dplot(res_our$u_mat[,1:3], cluster_labels,
                   bg_col_vec = col_vec2,
                   breaks = seq(0.5, 13.5, by = 1),
                   cluster_center = cluster_center,
                   center_col_vec = col_vec,
                   center_labels = 1:13,
                   curves = NA,
                   pch = 16, main = "eSVD embedding and trajectories\n(Curved Gaussian)",
                   xlab = "", ylab = "",
                   zlab = "Latent dimension 3",
                   theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                   xlim = bound_matrix[[kk]]$xlim,
                   ylim = bound_matrix[[kk]]$ylim,
                   zlim = bound_matrix[[kk]]$zlim)

  curves <- our_curves$curves
  col_vec_short <- color_func(0.9)[c(6,4)]
  for(i in 1:length(curves)){
    ord <- curves[[i]]$ord
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = "black", lwd = 2)
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = col_vec_short[i], lwd = 2)
  }
  graphics.off()
}

# 3D plots with tubes

our_tube_list <- lapply(1:length(our_curves$curves), function(x){
  s_mat <- our_curves$curves[[x]]$s[our_curves$curves[[x]]$ord,]
  construct_3d_tube(s_mat, radius = our_sd_val$sd_val)
})


for(kk in 1:nrow(angle_matrix)){
  print(kk)

  png(paste0("../figure/main/eSVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], "_tube.png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0,0,4,0))
  slingshot_3dplot(res_our$u_mat[,1:3], cluster_labels,
                   bg_col_vec = col_vec2,
                   breaks = seq(0.5, 13.5, by = 1),
                   cluster_center = cluster_center,
                   center_col_vec = col_vec,
                   center_labels = 1:13,
                   curves = NA,
                   pch = 16, main = "eSVD embedding with uncertainty tubes\n(Curved Gaussian)",
                   xlab = "", ylab = "",
                   zlab = "Latent dimension 3",
                   theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                   xlim = bound_matrix[[kk]]$xlim,
                   ylim = bound_matrix[[kk]]$ylim,
                   zlim = bound_matrix[[kk]]$zlim)

  curves <- our_curves$curves
  col_vec_short <- color_func(0.9)[c(6,4)]
  for(i in 1:length(curves)){
    ord <- curves[[i]]$ord
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = "black", lwd = 2)
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = col_vec_short[i], lwd = 2)

    col_mat <- our_tube_list[[i]]$z_mat

    plot3D::surf3D(our_tube_list[[i]]$x_mat,
                   our_tube_list[[i]]$y_mat,
                   our_tube_list[[i]]$z_mat, add = T,
                   colvar = col_mat,
                   col = colorRampPalette(c("white", col_vec_short[i]))(100),
                   breaks = seq(min(col_mat), max(col_mat), length.out = 101),
                   colkey = F)
  }
  graphics.off()
}

#######################################

# now for the svd embedding

svd_res <- svd(dat_impute)
svd_u <- svd_res$u[,1:p] %*% diag(sqrt(svd_res$d[1:p]))

cluster_center <- .compute_cluster_center(svd_u[,1:3], .construct_cluster_matrix(cluster_labels))

num_order_vec <- c(5, rep(3,2), rep(1,6), rep(2,2),  rep(5,2))
col_vec <- color_func(1)[num_order_vec]
col_vec2 <- color_func(0.1)[num_order_vec]
col_name <- c("orange", rep("bluish green", 2), rep("yellow", 6), rep("skyblue", 2), rep("orange", 2))
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                       level = sapply(1:13, function(y){which(sapply(cluster_group_list, function(x){y %in% x}))}),
                       col_name = col_name,
                       col_code = col_vec)
col_info

#angle_matrix
angle_matrix = matrix(c(225,225, 90,135, 270,225), byrow = T, ncol = 2)
bound_matrix = list(list(xlim = c(-7, 0), ylim = c(-5,4), zlim = c(-2,10)),
                    list(xlim = c(-7, 0), ylim = c(-5,4), zlim = c(-2,10)),
                    list(xlim = c(-7, 0), ylim = c(-5,4), zlim = c(-2,10)))

for(kk in 1:nrow(angle_matrix)){
  print(kk)

  # svd 3D plots
  png(paste0("../figure/main/SVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(1,1,4,1))
  slingshot_3dplot(svd_u[,1:3], cluster_labels,
                   bg_col_vec = col_vec2,
                   breaks = seq(0.5, 13.5, by = 1),
                   cluster_center = cluster_center,
                   center_col_vec = col_vec,
                   center_labels = 1:13,
                   curves = NA,
                   pch = 16, main = "SVD embedding and trajectories\n(Constant-variance Gaussian)",
                   xlab = "", ylab = "",
                   zlab = "Latent dimension 3",
                   theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                   xlim = bound_matrix[[kk]]$xlim,
                   ylim = bound_matrix[[kk]]$ylim,
                   zlim = bound_matrix[[kk]]$zlim)

  curves <- naive_curves$curves
  col_vec_short <- color_func(0.9)[c(5,4)]
  for(i in 1:length(curves)){
    ord <- curves[[i]]$ord
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = "black", lwd = 2)
  }
  graphics.off()
}

# svd 3D plots with tubes

naive_tube_list <- lapply(1:length(naive_curves$curves), function(x){
  s_mat <- naive_curves$curves[[x]]$s[naive_curves$curves[[x]]$ord,]
  construct_3d_tube(s_mat, radius = naive_sd_val$sd_val)
})

for(kk in 1:nrow(angle_matrix)){
  print(kk)
  col_vec_short <- color_func(0.9)[c(4)]

  png(paste0("../figure/main/SVD_theta", angle_matrix[kk,1], "_phi", angle_matrix[kk,2], "_tube.png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(1,1,4,1))
  slingshot_3dplot(svd_u[,1:3], cluster_labels,
                   bg_col_vec = col_vec2,
                   breaks = seq(0.5, 13.5, by = 1),
                   cluster_center = cluster_center,
                   center_col_vec = col_vec,
                   center_labels = 1:13,
                   curves = NA,
                   pch = 16, main = "SVD embedding and trajectories\n(Constant-variance Gaussian)",
                   xlab = "", ylab = "",
                   zlab = "Latent dimension 3",
                   theta = angle_matrix[kk,1], phi = angle_matrix[kk,2],
                   xlim = bound_matrix[[kk]]$xlim,
                   ylim = bound_matrix[[kk]]$ylim,
                   zlim = bound_matrix[[kk]]$zlim)

  curves <- naive_curves$curves
  col_vec_short <- color_func(0.9)[c(4)]
  for(i in 1:length(curves)){
    ord <- curves[[i]]$ord
    plot3D::lines3D(x = curves[[i]]$s[ord, 1],
                    y = curves[[i]]$s[ord, 2],
                    z = curves[[i]]$s[ord, 3],
                    add = T, colkey = F, col = "black", lwd = 2)
  }

  i <- 3
  col_mat <- naive_tube_list[[i]]$z_mat

  plot3D::surf3D(naive_tube_list[[i]]$x_mat,
                 naive_tube_list[[i]]$y_mat,
                 naive_tube_list[[i]]$z_mat, add = T,
                 colvar = col_mat,
                 col = colorRampPalette(c("white", col_vec_short[1]))(100),
                 breaks = seq(min(col_mat), max(col_mat), length.out = 101),
                 colkey = F)
  graphics.off()

}

