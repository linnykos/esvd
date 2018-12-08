load("../results/step4_factorization.RData")

.b_estimate <- function(mat, cluster_labels){
  uniq_val <- sort(unique(cluster_labels))
  mean(sapply(uniq_val, function(x){
    idx <- which(cluster_labels == x)
    stats::princomp(mat[idx,])$sdev[1]
  }))/4
}

k_select <- 3
u_mat <- res_our$u_mat[,1:k_select]

cluster_labels <- singlecell::dbscan(u_mat, neighbor_count = 10, upper_cutoff = 14,
                     size_cutoff = 20)

length(which(is.na(cluster_labels)))/length(cluster_labels)

upscale_vec <- max(table(cluster_labels))/table(cluster_labels)
# curves_list <- vector("list", length = length(table(cluster_labels)))
# for(i in 1:length(table(cluster_labels))){
#   curves_list[[i]] <- singlecell::slingshot(u_mat, cluster_labels, starting_cluster = i, b = 0.5, shrink = 1)
# }
# sapply(curves_list, function(x){length(x$lineages)})
# curves <- curves_list[[which.min(sapply(curves_list, function(x){length(x$lineages)}))]]

curves <- singlecell::slingshot(u_mat, cluster_labels, starting_cluster = 4, b = 0.5, shrink = 1,
                                upscale_vec = upscale_vec)

sapply(curves$lineages, function(x){
  length(which(cluster_labels %in% x))/length(cluster_labels)
})


###################

# SVD embedding
u_mat_svd <- res_naive$u[,1:k_select] %*% sqrt(diag(res_naive$d[1:k_select]))

cluster_labels_svd <- singlecell::dbscan(u_mat_svd, neighbor_count = 10, upper_cutoff = 14,
                                     size_cutoff = 45)
table(cluster_labels_svd)

col_vec <- cluster_labels_svd; col_vec[is.na(cluster_labels_svd)] <- 0; col_vec <- col_vec+1
plot(u_mat_svd[,1], u_mat_svd[,2], asp = T, pch = 16, col = col_vec)

# curves_list_svd <- vector("list", length = length(table(cluster_labels_svd)))
# for(i in 1:length(table(cluster_labels_svd))){
#   curves_list_svd[[i]] <-  .get_lineages(u_mat_svd, cluster_labels_svd, starting_cluster = i,
#                                          knn = NA, remove_outlier = F)
# }
#
# sapply(curves_list_svd, function(x){length(x)})

b_est <- .b_estimate(u_mat_svd, cluster_labels_svd)
curves_svd <- singlecell::slingshot(u_mat_svd, cluster_labels_svd, starting_cluster = 3,
                                    b = b_est, shrink = 1)

print(paste0(Sys.time(), ": Finished clustering"))
save.image(paste0("../results/step5_clustering", suffix, ".RData"))

##################
#
col_vec <- cluster_labels
col_vec[is.na(col_vec)] <- rgb(0,0,0,0.1)
plot(u_mat[,1], u_mat[,2], asp = T, pch = 16, col = col_vec)
plot(u_mat[,1], u_mat[,3], asp = T, pch = 16, col = col_vec)
plot(u_mat[,2], u_mat[,3], asp = T, pch = 16, col = col_vec)
c(length(is.na(cluster_labels)), table(cluster_labels))

# k <- 1
# par(mfrow = c(1,3))
# idx <- which(cluster_labels %in% curves$lineages[[k]])
# combn_mat <- utils::combn(3, 2)
# mid_vec <- apply(u_mat, 2, function(x){mean(range(x))})[1:3]
# rg <- max(apply(u_mat, 2, function(x){diff(range(x))})[1:3])
# lim_list <- lapply(1:3, function(x){mid_vec[x]+c(-1,1)*rg/2})
# col_vec <- cluster_labels
# col_vec[which(is.na(col_vec))] <- rgb(0,0,0)
# for(i in 1:ncol(combn_mat)){
#   idx1 <- combn_mat[1,i]; idx2 <- combn_mat[2,i]
#   plot(u_mat[-idx,idx1], u_mat[-idx,idx2], pch = 16, col = rgb(0.85,0.85,0.85),
#        asp = T, cex = 1.3, xlim = lim_list[[idx1]], ylim = lim_list[[idx2]],
#        xaxt = "n", yaxt = "n")
#   points(u_mat[idx,idx1], u_mat[idx,idx2], pch = 16, col = "white",
#          cex = 1.3)
#   points(u_mat[idx,idx1], u_mat[idx,idx2], pch = 16, col = col_vec[idx],
#          cex = 1)
# }
#
#
# ##
#
# for(i in 1:ncol(combn_mat)){
#   idx1 <- combn_mat[1,i]; idx2 <- combn_mat[2,i]
#   plot(u_mat[,idx1], u_mat[,idx2], pch = 16, col = rgb(0.85,0.85,0.85),
#        asp = T, cex = 1.3, xlim = lim_list[[idx1]], ylim = lim_list[[idx2]],
#        xlab = paste0("Latent dimension ", idx1),
#        ylab = paste0("Latent dimension ", idx2))
#
#   for(k in 1:length(curves$curves)){
#     ord <- curves$curves[[k]]$ord
#     lines(curves$curves[[k]]$s[ord, idx1], curves$curves[[k]]$s[ord, idx2], lwd = 3.5,
#           col = "white")
#     lines(curves$curves[[k]]$s[ord, idx1], curves$curves[[k]]$s[ord, idx2], lwd = 3,
#           col = "black")
#   }
#
#   cluster_mat <- .construct_cluster_matrix(cluster_labels)
#   centers <- .compute_cluster_center(u_mat, cluster_mat)
#   points(centers[,idx1], centers[,idx2], col = "white", pch = 16, cex = 2.25)
#   points(centers[,idx1], centers[,idx2], col = 1:11, pch = 16, cex = 2)
# }


