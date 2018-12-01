load("../results/step4_clustering.RData")

u_mat <- res$u_mat[,c(2,1,3)]
col_vec <- cluster_labels
col_vec[which(is.na(col_vec))] <- rgb(0,0,0)

col_template <- c(rgb(144,113,38, maxColorValue=255), #brown
                  rgb(0,0,0),
                  rgb(149,219,144, maxColorValue=255), #green
                  rgb(227,73,86, maxColorValue=255), #red
                  rgb(100,100,200, maxColorValue=255), #purple
                  rgb(40,225,201, maxColorValue=255), #turqouise
                  rgb(100,140,252, maxColorValue=255), #blue
                  rgb(255,152,41, maxColorValue=255), #orange
                  rgb(238,204,17, maxColorValue=255) #goldenrod
)

for(i in 1:max(cluster_labels, na.rm = T)){
  col_vec[col_vec == i] <- col_template[i]
}

########################

# plot cell trajectories
combn_mat <- utils::combn(3, 2)
mid_vec <- apply(u_mat, 2, function(x){mean(range(x))})[1:3]
rg <- max(apply(u_mat, 2, function(x){diff(range(x))})[1:3])
lim_list <- lapply(1:3, function(x){mid_vec[x]+c(-1,1)*rg/2})

png("../figure/main/latent_cluster.png", height = 1500, width = 2400, res = 300, units = "px")
par(mfcol = c(3,length(curves$lineages)), mar = c(0.5,0.5,0.5,0.5))
for(k in c(2,1,4,3,5)){
  idx <- which(cluster_labels %in% curves$lineages[[k]])

  for(i in 1:ncol(combn_mat)){
    idx1 <- combn_mat[1,i]; idx2 <- combn_mat[2,i]
    plot(u_mat[-idx,idx1], u_mat[-idx,idx2], pch = 16, col = rgb(0.85,0.85,0.85),
         asp = T, cex = 1.3, xlim = lim_list[[idx1]], ylim = lim_list[[idx2]],
         xaxt = "n", yaxt = "n")
    points(u_mat[idx,idx1], u_mat[idx,idx2], pch = 16, col = "white",
           cex = 1.3)
    points(u_mat[idx,idx1], u_mat[idx,idx2], pch = 16, col = col_vec[idx],
           cex = 1)
  }
}
graphics.off()

########

# trajectories
combn_mat <- utils::combn(3, 2)
mid_vec <- apply(u_mat, 2, function(x){mean(range(x))})[1:3]
rg <- max(apply(u_mat, 2, function(x){diff(range(x))})[1:3])
lim_list <- lapply(1:3, function(x){mid_vec[x]+c(-1,1)*rg/2})

png("../figure/main/trajectory.png", height = 800, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3), mar = c(4,4,0.5,0.5))
for(i in 1:ncol(combn_mat)){
  idx1 <- combn_mat[1,i]; idx2 <- combn_mat[2,i]
  plot(u_mat[,idx1], u_mat[,idx2], pch = 16, col = rgb(0.85,0.85,0.85),
       asp = T, cex = 1.3, xlim = lim_list[[idx1]], ylim = lim_list[[idx2]],
       xlab = paste0("Latent dimension ", idx1),
       ylab = paste0("Latent dimension ", idx2))

  for(k in 1:length(curves$curves)){
    ord <- curves$curves[[k]]$ord
    lines(curves$curves[[k]]$s[ord, idx1], curves$curves[[k]]$s[ord, idx2], lwd = 3.5,
          col = "white")
    lines(curves$curves[[k]]$s[ord, idx1], curves$curves[[k]]$s[ord, idx2], lwd = 3,
          col = "black")
  }

  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  centers <- .compute_cluster_center(u_mat, cluster_mat)
  points(centers[,idx1], centers[,idx2], col = "white", pch = 16, cex = 2.25)
  points(centers[,idx1], centers[,idx2], col = col_template, pch = 16, cex = 2)
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
