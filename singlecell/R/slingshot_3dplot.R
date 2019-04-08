slingshot_3dplot <- function(dat, cluster_labels, col_vec, breaks, curves, ...){
  stopifnot(length(col_vec) == length(breaks) + 1)

  cluster_center <- .compute_cluster_center(dat, .construct_cluster_matrix(cluster_labels))

  plot3D::scatter3D(x = dat[,1], y = dat[,2], z = dat[,3],
                    surface = FALSE, colvar = cluster_labels,
                    breaks = breaks, col = col_vec, ...)
  plot3D::points3D(cluster_center[,1], cluster_center[,2], cluster_center[,3],
                   colvar = 1:max(cluster_labels),
                   breaks = breaks, add = T, col = col_Vec, colkey = F, ...)

  for(k in 1:length(curves$curves)){
    ord <- curves$curves[[k]]$ord
    plot3D::lines3D(x = curves$curves[[k]]$s[ord, 1],
                    y = curves$curves[[k]]$s[ord, 2],
                    z = curves$curves[[k]]$s[ord, 3],
                    add = T, lwd = 2, colkey = F, col = "black")
  }

  invisible()
}

slingshot_construct_tube <- function(curves, sd_vec){

}

##############

slingshot_construct_tube_singleton <- function(curve, sd_val){

}
