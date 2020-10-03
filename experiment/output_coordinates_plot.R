rm(list=ls())
embedding <- read.csv("../results/embedding_coordinates.csv")
curve1 <- read.csv("../results/curve1_coordinates.csv")
curve2 <- read.csv("../results/curve2_coordinates.csv")

head(embedding)
head(curve1)
head(curve2)
curve_list <- list(curve1, curve2)

###############################

color_func <- function(alpha1, alpha2){
  c(grDevices::rgb(240/255, 228/255, 66/255, alpha1), #yellow
    grDevices::rgb(86/255, 180/255, 233/255, alpha1), #skyblue
    grDevices::rgb(0/255, 158/255, 115/255, alpha1), #bluish green
    grDevices::rgb(0/255, 114/255, 178/255,alpha1), #blue
    grDevices::rgb(230/255, 159/255, 0/255,alpha1), #orange
    grDevices::rgb(100/255, 100/255, 100/255, alpha2)) #gray
}

cluster_sub_vec <- sort(as.character(unique(embedding$cell_subtype)))
# hand-coded colors
color_vec_bg <- color_func(0.2, 0.1)[c(5, rep(3,2), c(6,1,1,4,4,4), rep(2,2), rep(5,2))]
color_vec_center <- color_func(1, 1)[c(5, rep(3,2), c(6,1,1,4,4,4), rep(2,2), rep(5,2))]
center_mat <- t(sapply(cluster_sub_vec, function(x){
  idx <- which(embedding$cell_subtype == x)
  colMeans(embedding[idx,-c(1:2)])
}))
color_table <- data.frame(cell_subtype = cluster_sub_vec, color_bg = color_vec_bg,
                          color_center = color_vec_center,
                          center = center_mat)

theta <- 270; phi <- 225
spacing <- 4.5; tmp <- spacing*c(-.5,.5)
xlim <- -3+tmp; ylim <- 0.75+tmp; zlim <- 1+tmp

###################################

breaks <- seq(0.5, nrow(color_table)+1, by = 1)
par(mar = c(0,0,4,0))

# plot 3d points
plot3D::scatter3D(x = embedding$coord.1, y = embedding$coord.2, z = embedding$coord.3,
                  surface = FALSE,
                  colvar = as.numeric(sapply(embedding$cell_subtype, function(x){which(color_table$cell_subtype == x)[1]})),
                  breaks = breaks, col = color_table$color_bg, colkey = F, theta = theta, phi = phi,
                  cex = 0.4, pch = 16, xlim = xlim, ylim = ylim, zlim = zlim,
                  main = "eSVD embedding and trajectories:\nCurved Gaussian",
                  xlab = "Latent dimension 1", ylab = "Latent dimension 2",
                  zlab = "Latent dimension 3")

# plot centers
plot3D::points3D(x = color_table$center.coord.1, y = color_table$center.coord.2, z = color_table$center.coord.3,
                 colvar =  1:nrow(color_table),
                 breaks = breaks, add = T, col = rep("black", 13),
                 cex = 1.5*2, colkey = F, pch = 16)
plot3D::points3D(x = color_table$center.coord.1, y = color_table$center.coord.2, z = color_table$center.coord.3,
                 colvar =  1:nrow(color_table),
                 breaks = breaks, add = T, col = color_table$color_center,
                 cex = 2, colkey = F, pch = 16)

# plot trajectories
col_vec_short <- color_func(1,1)[c(1,4)]
for(i in 1:length(curve_list)){
  plot3D::lines3D(x = curve_list[[i]][, 1],
                  y = curve_list[[i]][, 2],
                  z = curve_list[[i]][, 3],
                  add = T, colkey = F, col = "black", lwd = 12)
  plot3D::lines3D(x = curve_list[[i]][, 1],
                  y = curve_list[[i]][, 2],
                  z = curve_list[[i]][, 3],
                  add = T, colkey = F, col = col_vec_short[i], lwd = 12)

}

