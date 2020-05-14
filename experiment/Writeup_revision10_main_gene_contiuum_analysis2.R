rm(list=ls())
load("../results/step5_trajectory_original.RData")

upscale_factor <- 1
reduction_percentage <- 0.2

p <- 3
set.seed(10)
esvd_curves <- eSVD::slingshot(esvd_embedding$u_mat[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                               cluster_group_list = cluster_group_list,
                               verbose = T, upscale_factor = upscale_factor,
                               reduction_percentage = reduction_percentage,
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

slingshot_res <- esvd_curves
pseudotime_df_list <- .extract_pseudotimes(slingshot_res)
shared_df <- .compile_common_cells(pseudotime_df_list)

table(shared_df$consensus, cluster_labels[shared_df$cell_idx])

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

plot(NA, xlim = c(0, max(shared_df$pseudotime)), ylim = c(0,7))
for(i in 1:6){
  tmp <- shared_df$pseudotime[row_idx[shared_df$cell_idx] == i]
  den_res <- density(tmp)
  y_vec <- (c(0, den_res$y, 0 , 0))
  y_vec <- y_vec/1.5
  polygon(x = c(den_res$x[1], den_res$x, den_res$x[length(den_res$x)], den_res$x[1]),
          y = y_vec + i,
          col = col_vec_svd[which(floor(order_vec_svd) == i)][1])
  # points(tmp, rep(i, length(tmp)), pch = 1, col = col_vec2_svd[which(floor(order_vec_svd) == i)][1])
}

plot(shared_df$pseudotime, col = col_vec2_svd[cluster_labels[shared_df$cell_idx]], pch = 16)
shared_df2 <- shared_df[which(shared_df$consensus),]
plot(shared_df2$pseudotime, col = col_vec2_svd[cluster_labels[shared_df2$cell_idx]], pch = 16)

