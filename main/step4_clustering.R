load("../results/step3_factorization_logged.RData")

u_mat <- res$u_mat[,1:3]
cluster_labels <- dbscan(u_mat, neighbor_count = 10, upper_cutoff = 14,
                     size_cutoff = 19)
table(cluster_labels)

col_vec <- cluster_labels
col_vec[which(is.na(col_vec))] <- rgb(0,0,0,0.1)
plot(u_mat[,1], u_mat[,2], pch = 16, col = col_vec, asp = T)
plot(u_mat[,1], u_mat[,3], pch = 16, col = col_vec, asp = T)
plot(u_mat[,2], u_mat[,3], pch = 16, col = col_vec, asp = T)

# plot3D::scatter3D(u_mat[,1], u_mat[,2], u_mat[,3], col = col_vec,
#                   pch = 16)

############

# determine the lineage
lineages <- .get_lineages(u_mat, cluster_labels, starting_cluster = 1, knn = NA,
              remove_outlier = T, percentage = 0.05)

