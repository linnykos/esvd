rm(list=ls())
load("../results/step4_factorization_spca.RData")
dat <- res_our$u_mat
starting_cluster <- 1

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

dist_mat <- .compute_cluster_distances(dat, cluster_labels)
starting_cluster <- cluster_group_list[[1]][1]

lineages <- .get_lineages(dat, cluster_labels, starting_cluster = starting_cluster,
                          cluster_group_list = cluster_group_list)

intersect_1 <- setdiff(lineages[[1]], lineages[[2]])
intersect_2 <- setdiff(lineages[[2]], lineages[[1]])
factor_vec <- rep(1, length(cell_type_vec))
factor_vec[which(cluster_labels %in% intersect_1)] <- 2
factor_vec[which(cluster_labels %in% intersect_2)] <- 3

# let's try a 3d plot
car::scatter3d(x = res_our$u_mat[,1], y = res_our$u_mat[,2], z = res_our$u_mat[,3],
               surface=FALSE, groups = as.factor(factor_vec))


##################
shrink = 1
thresh = 0.001
max_iter = 15
b = 1
upscale_vec = NA
curves <- .get_curves(dat, cluster_labels, lineages, shrink = shrink,
                      thresh = thresh, max_iter = max_iter, b = b, upscale_vec = upscale_vec)

idx1 <- 2; idx2 <- 3
plot(res_our$u_mat[,idx1], res_our$u_mat[,idx2], asp = T, col = "gray", pch = 16)

for(k in 1:2){
  ord <- curves[[k]]$ord
  lines(curves[[k]]$s[ord, idx1], curves[[k]]$s[ord, idx2], lwd = 3.5,
        col = "white")
  lines(curves[[k]]$s[ord, idx1], curves[[k]]$s[ord, idx2], lwd = 3,
        col = "black")
}
