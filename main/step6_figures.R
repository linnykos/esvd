load("../results/step5_clustering.RData")

col_vec <- cluster_labels
col_vec[which(is.na(col_vec))] <- rgb(0,0,0)

col_template <- c(rgb(238,204,17, maxColorValue=255), #goldenrod
                  rgb(149,219,144, maxColorValue=255), #green
                  rgb(0,0,0),
                  rgb(227,73,86, maxColorValue=255), #red
                  rgb(100,100,200, maxColorValue=255), #purple
                  rgb(40,225,201, maxColorValue=255), #turqouise
                  rgb(100,140,252, maxColorValue=255), #blue
                  rgb(255,152,41, maxColorValue=255), #orange
                  rgb(144,113,38, maxColorValue=255) #brown
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

cluster_mat <- .construct_cluster_matrix(cluster_labels)
centers <- .compute_cluster_center(u_mat, cluster_mat)

png("../figure/main/latent_cluster.png", height = 1200, width = 1900, res = 300, units = "px")
par(mfcol = c(3,length(curves$lineages)), mar = c(0.5,0.5,0.5,0.5))
for(k in c(2,1,5,4,3)){
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

    # ord <- curves$curves[[k]]$ord
    # lines(curves$curves[[k]]$s[ord, idx1], curves$curves[[k]]$s[ord, idx2], lwd = 6,
    #       col = "white")
    # lines(curves$curves[[k]]$s[ord, idx1], curves$curves[[k]]$s[ord, idx2], lwd = 4,
    #       col = "black")
    #
    # points(centers[curves$lineages[[k]], idx1], centers[curves$lineages[[k]], idx2],
    #        col = "white", pch = 16, cex = 2.5)
    # points(centers[curves$lineages[[k]], idx1], centers[curves$lineages[[k]], idx2],
    #        col = col_template[curves$lineages[[k]]], pch = 16, cex = 2)
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
for(i in 1:3){
  idx1 <- combn_mat[1,i]; idx2 <- combn_mat[2,i]
  plot(u_mat[,idx1], u_mat[,idx2], pch = 16, col = rgb(0.85,0.85,0.85,1),
       asp = T, cex = 1, xlim = lim_list[[idx1]], ylim = lim_list[[idx2]],
       xlab = paste0("Latent dimension ", idx1),
       ylab = paste0("Latent dimension ", idx2))

  for(k in 1:length(curves$lineages)){
    ord <- curves$curves[[k]]$ord
    lines(curves$curves[[k]]$s[ord, idx1], curves$curves[[k]]$s[ord, idx2], lwd = 3.5,
          col = "white")
    lines(curves$curves[[k]]$s[ord, idx1], curves$curves[[k]]$s[ord, idx2], lwd = 3,
          col = "black")
  }

  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  centers <- .compute_cluster_center(u_mat, cluster_mat)
  points(centers[,idx1], centers[,idx2], col = "white", pch = 16, cex = 1.75)
  points(centers[,idx1], centers[,idx2], col = col_template, pch = 16, cex = 1.5)
}
graphics.off()

#############

# plot the diagnostic
idx <- which(scalar_vec == scalar_val)
pred_mat <- 1/(res_list[[idx]]$u_mat %*% t(res_list[[idx]]$v_mat))
pred_mat <- t(sapply(1:nrow(pred_mat), function(x){
  pred_mat[x,] * extra_weight[x]
}))
xlim <- range(c(pred_mat[missing_idx], pred_naive[missing_idx], dat_impute[missing_idx]))
ylim <- xlim

diag_vec <- c(0,1)
mat <- cbind(dat_impute[missing_idx], pred_mat[missing_idx])
our_pca <- stats::princomp(mat)
our_angle <- as.numeric(acos(diag_vec %*% our_pca$loadings[,1]))
mat <- cbind(dat_impute[missing_idx], pred_naive[missing_idx])
naive_pca <- stats::princomp(mat)
naive_angle <- as.numeric(acos(diag_vec %*% naive_pca$loadings[,1]))

rad <- 500

png("../figure/main/diagnostic.png", height = 1200, width = 2200, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,4,4,0.5))
plot(NA, asp = T, xlim = xlim, ylim = ylim, xlab = "Predicted mean",
     ylab = "Observed value", main = "Gaussian model with\nfixed variance")
lines(c(-1e6,1e6), rep(0,2), col = "red", lwd = 2)
lines(rep(0,2), c(-1e6,1e6), col = "red", lwd = 2)
points(pred_naive[missing_idx], dat_impute[missing_idx], pch = 16, col = rgb(0,0,0,0.2))
lines(c(-1e6,1e6), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, 1e6), c(0, 1e6*naive_pca$loadings[1,1]/naive_pca$loadings[2,1]), col = "blue", lwd = 2, lty = 2)
theta_seq <- seq(0, naive_angle, length.out = 100)
lines(rad*cos(theta_seq), rad*sin(theta_seq), lwd = 1.5, lty = 2)
text(rad,200, labels = paste0(round(naive_angle*180/pi,1), " degrees"), pos = 4)

plot(NA, asp = T, xlim = xlim, ylim = ylim, xlab = "Predicted mean",
     ylab = "Observed value", main = "Curved Gaussian model")
lines(c(-1e6,1e6), rep(0,2), col = "red", lwd = 2)
lines(rep(0,2), c(-1e6,1e6), col = "red", lwd = 2)
points(pred_mat[missing_idx], dat_impute[missing_idx], pch = 16, col = rgb(0,0,0,0.2))
lines(c(-1e6,1e6), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, 1e6), c(0, 1e6*our_pca$loadings[1,1]/our_pca$loadings[2,1]), col = "blue", lwd = 2, lty = 2)
theta_seq <- seq(0, our_angle, length.out = 100)
lines(rad*cos(theta_seq), rad*sin(theta_seq), lwd = 1.5, lty = 2)
text(rad,200, labels = paste0(round(our_angle*180/pi,1), " degrees"), pos = 4)
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
