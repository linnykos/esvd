rm(list=ls())
load("../results/step5_trajectory_original.RData")

p <- 3
set.seed(10)
esvd_curves <- eSVD::slingshot(esvd_embedding$u_mat[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                               cluster_group_list = cluster_group_list, shrink = 2,
                               verbose = T, upscale_factor = 1, stretch = 9999, max_iter = 3,
                               squared = T)

# checking of all points are accounted for
idx <- which(cluster_labels %in% esvd_curves$lineages[[1]])
## first extract which indices are in the first lineage
idx2 <- which(esvd_curves$curves[[1]]$W != 0)
all(idx %in% esvd_curves$idx[idx2])

## same thing for second lineage
idx <- which(cluster_labels %in% esvd_curves$lineages[[2]])
idx2 <- which(esvd_curves$curves[[2]]$W != 0)
all(idx %in% esvd_curves$idx[idx2])

###################################

pseudotime_df <- construct_pseudotime_trajectory_matrix(esvd_curves, cluster_labels)

pseudotime_df2 <- pseudotime_df[-intersect(which(is.na(pseudotime_df$consensus)), which(pseudotime_df$pseudotime <= 3)),]
pseudotime_df2 <- pseudotime_df2[-which(!pseudotime_df2$consensus),]
pseudotime_df2 <- pseudotime_df2[-intersect(which(pseudotime_df2$dist_to_curve >= 0.1),
                                            which(pseudotime_df2$status == "1")),]
zz_idx <- pseudotime_df2$cell_idx

# verify that cell type is roughly aligned with trajectory
color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(100/255, 100/255, 100/255, alpha)) #gray
}
num_order_vec_svd <- c(5, rep(3,2), c(1,1,1,1,1,1), rep(2,2),  rep(5,2))
col_vec_svd <- color_func(1)[num_order_vec_svd]
col_vec2_svd <- color_func(0.5)[num_order_vec_svd]
order_vec_svd <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
row_idx <- floor(order_vec_svd)[cluster_labels]

plot(NA, xlim = c(0, max(pseudotime_df$pseudotime)), ylim = c(0,7))
for(i in 1:6){
  tmp <- pseudotime_df$pseudotime[row_idx[pseudotime_df$cell_idx] == i]
  den_res <- density(tmp)
  y_vec <- (c(0, den_res$y, 0 , 0))
  y_vec <- y_vec/1.5
  polygon(x = c(den_res$x[1], den_res$x, den_res$x[length(den_res$x)], den_res$x[1]),
          y = y_vec + i,
          col = col_vec_svd[which(floor(order_vec_svd) == i)][1])
}

plot(NA, xlim = c(0, max(pseudotime_df2$pseudotime)), ylim = c(0,7))
for(i in 1:6){
  tmp <- pseudotime_df2$pseudotime[row_idx[pseudotime_df2$cell_idx] == i]
  den_res <- density(tmp)
  y_vec <- (c(0, den_res$y, 0 , 0))
  y_vec <- y_vec/1.5
  polygon(x = c(den_res$x[1], den_res$x, den_res$x[length(den_res$x)], den_res$x[1]),
          y = y_vec + i,
          col = col_vec_svd[which(floor(order_vec_svd) == i)][1])
}

## try another plot
plot(pseudotime_df$pseudotime, col = col_vec2_svd[cluster_labels[pseudotime_df$cell_idx]], pch = 16)
plot(pseudotime_df2$pseudotime, col = col_vec2_svd[cluster_labels[pseudotime_df2$cell_idx]], pch = 16)
# table(pseudotime_df2$cluster_labels[intersect(which(is.na(pseudotime_df2$consensus)), which(pseudotime_df2$pseudotime <= 3))])

## try another plot
col_palatte <- colorRampPalette(c("azure","darkviolet"))(100)
max_val <- max(pseudotime_df$pseudotime)
min_val <- min(pseudotime_df$pseudotime)
col_vec <-  col_palatte[pmin(pmax(round((pseudotime_df$pseudotime-min_val)/(max_val - min_val) * 100), 1), 100)]

cluster_center_esvd <- .compute_cluster_center(esvd_embedding$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))

num_order_vec_esvd <- c(5, rep(3,2), c(1,1,1,1,6,4), rep(2,2),  rep(5,2))
col_vec_esvd <- color_func(1)[num_order_vec_esvd]
col_vec2_esvd <- color_func(0.5)[num_order_vec_esvd]
col_name_esvd <- c("orange", rep("bluish green", 2), c("yellow", "yellow", "yellow", "yellow", "gray", "blue"), rep("skyblue", 2), rep("orange", 2))
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
plotting_order_esvd <- c(3,4,5,2,6,1)

combn_mat <- combn(3,2)
par(mfrow = c(1,3))
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding and trajectories\n(Curved Gaussian)")

  points(x = esvd_embedding$u_mat[zz_idx,i], y = esvd_embedding$u_mat[zz_idx,j], pch = 16, col = col_vec[zz_idx])

  for(ll in 1:nrow(cluster_center_esvd)){
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
  }

  curves <- esvd_curves$curves
  for(ll in 1:length(curves)) {
    ord <- curves[[ll]]$ord
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 8)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black", lwd = 5)
  }
}

# only one trajectory
traj_clust <- c(4,5,6,7)
# traj_clust <- 9
combn_mat <- combn(3,2)
par(mfrow = c(1,3))
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding and trajectories\n(Curved Gaussian)")

  points(x = esvd_embedding$u_mat[zz_idx[which(pseudotime_df2$cluster_labels %in% traj_clust)],i],
         y = esvd_embedding$u_mat[zz_idx[which(pseudotime_df2$cluster_labels %in% traj_clust)],j],
         pch = 16, col = col_vec[zz_idx[which(pseudotime_df2$cluster_labels %in% traj_clust)]])

  for(ll in 1:nrow(cluster_center_esvd)){
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
  }

  curves <- esvd_curves$curves
  for(ll in 1:length(curves)) {
    ord <- curves[[ll]]$ord
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 8)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black", lwd = 5)
  }
}

# full plot

combn_mat <- combn(3,2)
par(mfrow = c(1,3))
for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding and trajectories\n(Curved Gaussian)")

  points(x = esvd_embedding$u_mat[,i], y = esvd_embedding$u_mat[,j], pch = 16, col = col_vec)

  for(ll in 1:nrow(cluster_center_esvd)){
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[ll])
  }

  curves <- esvd_curves$curves
  for(ll in 1:length(curves)) {
    ord <- curves[[ll]]$ord
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 8)
    lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black", lwd = 5)
  }
}

