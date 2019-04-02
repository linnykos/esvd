rm(list=ls())
load("../experiment/Week36_marques_slingshot.RData")

col_func <- function(alpha){
  col_vec <- numeric(length(unique(cluster_labels)))
  for(i in 1:3){
    col_vec[cluster_group_list[[i]]] <- rgb(238/255,204/255,17/255,alpha) #goldenrod
  }
  col_vec[cluster_group_list[[4]]] <- rgb(129/255,199/255,124/255,alpha) #green
  col_vec[cluster_group_list[[5]]] <- rgb(227/255,73/255,86/255,alpha) #red
  col_vec[cluster_group_list[[6]]] <- rgb(100/255,140/255,252/255,alpha) #blue

  col_vec[c(1,2,4,10,12)]
}


plot_3d_func <- function(dat, curves, name, theta = 60, phi = 20){
  cluster_center <- .compute_cluster_center(dat, .construct_cluster_matrix(cluster_labels))

  png(paste0("../figure/experiment/Week36_3dplots/Week36_marques_", name,
             "_lineage_theta", theta, "_phi", phi, ".png"),
      height = 2000, width = 2000, res = 300, units = "px")
  par(mar = c(0.5, 0.5, 4, 0.5))
  plot3D::scatter3D(x = dat[,1], y = dat[,2], z = dat[,3],
                    surface = FALSE, colvar = cluster_labels,
                    breaks = c(0, 1.5, 3.5, 9.5, 11.5, 14),
                    col = col_func(0.1), pch = 16, cex = 0.4, colkey = F,
                    theta = theta, phi = phi, main = paste0(name, " lineages"),
                    xlab = "Latent dimension 1",
                    ylab = "Latent dimension 2",
                    zlab = "Latent dimension 3")
  plot3D::points3D(cluster_center[,1], cluster_center[,2], cluster_center[,3],
                   colvar = 1:13,
                   breaks = c(0, 1.5, 3.5, 9.5, 11.5, 14),
                   add = T, pch = 16, cex = 2, col = col_func(1), colkey = F)

  reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*.25
  for(k in 1:length(curves$curves)){
    ord <- curves$curves[[k]]$ord
    plot3D::lines3D(x = curves$curves[[k]]$s[ord, 1]*reduction_factor,
                    y = curves$curves[[k]]$s[ord, 2]*reduction_factor,
                    z = curves$curves[[k]]$s[ord, 3]*reduction_factor,
                    add = T, lwd = 2, colkey = F, col = "black")
  }
  graphics.off()
}

plot_3d_func(u_mat, our_curves, "Our")

paramMat <- as.matrix(expand.grid(seq(0, 360, length.out = 17),
                                  seq(0, 360, length.out = 17)))
sapply(1:nrow(paramMat), function(x){
  if(x %% round(nrow(paramMat)/10) == 0) cat('*')
  plot_3d_func(u_mat, our_curves, "Our", theta = paramMat[x,1], phi = paramMat[x,2])
})

plot_3d_func(naive, naive_curves, "Naive")

sapply(1:nrow(paramMat), function(x){
  if(x %% round(nrow(paramMat)/10) == 0) cat('*')
  plot_3d_func(naive, naive_curves, "Naive", theta = paramMat[x,1], phi = paramMat[x,2])
})

#
# # information from https://rdrr.io/cran/car/man/scatter3d.html
# car::scatter3d(x = u_mat[,1], y = u_mat[,2], z = u_mat[,3],
#                surface=FALSE, point.col = col_vec[cluster_labels], pch = 1)
#
# # how to get angles https://stackoverflow.com/questions/22257196/get-rgl-view-parameters
