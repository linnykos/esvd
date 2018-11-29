load("../results/step4_clustering.RData")

u_mat <- res$u_mat[,1:k]
col_vec <- cluster_labels
col_vec[which(is.na(col_vec))] <- rgb(0,0,0,0.1)

col_template <- c(rgb(0,0,0,0.5),
                  rgb(205,40,54, 120, maxColorValue=255), #red
                  rgb(180,200,255, 230, maxColorValue=255), #blue
                  rgb(100,100,200, 200, maxColorValue=255), #purple
                  rgb(149,219,144, 200, maxColorValue=255), #green
                  rgb(238,204,17, 200, maxColorValue=255), #goldenrod
                  rgb(40,225,201, 200, maxColorValue=255), #turqouise
                  rgb(255,152,41, 150, maxColorValue=255) #orange
)

for(i in 1:max(cluster_labels, na.rm = T)){
  col_vec[col_vec == i] <- col_template[i]
}

########################

# plot cell trajectories

png("../figure/main/trajectory.png", height = 800, width = 2200, res = 300, units = "px")
par(mfrow = c(1,3), mar = c(4,4,0.5,0.5))
combn_mat <- utils::combn(3, 2)
for(i in 1:ncol(combn_mat)){
  idx1 <- combn_mat[1,i]; idx2 <- combn_mat[2,i]
  plot(u_mat[,idx1], u_mat[,idx2], pch = 16, col = col_vec, asp = T,
       xlab = paste0("Latent dimension ", idx1),
       ylab = paste0("Latent dimension ", idx2))
  for(i in 1:length(curves$curves)){
    ord <- curves$curves[[i]]$ord
    lines(curves$curves[[i]]$s[ord, idx1], curves$curves[[i]]$s[ord, idx2], lwd = 3.5,
          col = "white")
    lines(curves$curves[[i]]$s[ord, idx1], curves$curves[[i]]$s[ord, idx2], lwd = 3,
          col = "black")
  }
}
graphics.off()

##########################

# plot estimated density
png("../figure/main/estimated_density.png", height = 800, width = 2200, res = 300, units = "px")
par(mfrow = c(1,3), mar = c(4,4,0.5,0.5))
combn_mat <- utils::combn(3, 2)
h_vec <- c(0.4, 0.6, 0.6)
for(i in 1:ncol(combn_mat)){
  idx1 <- combn_mat[1,i]; idx2 <- combn_mat[2,i]
  den <- MASS::kde2d(u_mat[,idx1], u_mat[,idx2], h = h_vec[i], n = 100)
  image(den, col = grDevices::heat.colors(100, alpha = 0.7),
        xlab = paste0("Latent dimension ", idx1),
        ylab = paste0("Latent dimension ", idx2), asp = T)
  points(u_mat[,idx1], u_mat[,idx2], col = rgb(0, 0,0, 0.3), pch = 16)
  contour(den, add = T, drawlabels = F, col = rgb(1,1,1,0.5), lwd = 2,
          nlevels = 10)
}
graphics.off()

#######
# 3D plot

# plot3D::scatter3D(u_mat[,1], u_mat[,2], u_mat[,3], col = rgb(0,0,0,0), pch = 16,
#                   theta = 40, phi = 15)
# for(i in 1:length(curves$curves)){
#   ord <- curves$curves[[i]]$ord
#   plot3D::lines3D(curves$curves[[i]]$s[ord, 1], curves$curves[[i]]$s[ord, 2],
#                   curves$curves[[i]]$s[ord, 3], col = "black", add = T)
# }

######
# plot raw data in their clusters

# idx <- which(is.na(cluster_labels))
# k <- max(cluster_labels, na.rm = T)
# for(i in k:1){
#   idx <- c(which(cluster_labels == i), idx)
# }
#
# .plot_singlecell(dat[idx,])
# line_idx <- sapply(1:k, function(x){
#   1-length(which(cluster_labels <= x))/length(cluster_labels)
# })
# for(i in 1:k){
#   lines(x = c(0,1), y = rep(line_idx[i], 2), lwd = 2, lty = 2)
# }
#
