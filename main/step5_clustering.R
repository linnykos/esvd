load("../results/step4_factorization.RData")

k_select <- 3
u_mat <- res$u_mat[,1:k_select]; u_mat <- u_mat[,c(2,1,3)]
cluster_labels <- singlecell::dbscan(u_mat, neighbor_count = 10, upper_cutoff = 14,
                     size_cutoff = 26)

upscale_vec <- max(table(cluster_labels))/table(cluster_labels)
curves <- singlecell::slingshot(u_mat, cluster_labels, starting_cluster = 4, b = 1, shrink = 1,
                    upscale_vec = NA)

print(paste0(Sys.time(), ": Finished clustering"))
save.image("../results/step5_clustering.RData")

##################

# plot(res$u_mat[,1], res$u_mat[,2], asp = T, pch = 16, col = cluster_labels)
# plot(res$u_mat[,1], res$u_mat[,3], asp = T, pch = 16, col = cluster_labels)
# plot(res$u_mat[,2], res$u_mat[,3], asp = T, pch = 16, col = cluster_labels)
# c(length(is.na(cluster_labels)), table(cluster_labels))


# k <- 5
# par(mfrow = c(1,3))
# idx <- which(cluster_labels %in% curves$lineages[[k]])
# combn_mat <- utils::combn(3, 2)
# mid_vec <- apply(res$u_mat, 2, function(x){mean(range(x))})[1:3]
# rg <- max(apply(res$u_mat, 2, function(x){diff(range(x))})[1:3])
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


