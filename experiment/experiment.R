rm(list=ls())
load("../results/step4_factorization_spca.RData")
zz <- svd(dat)
naive <- zz$u[,1:3]%*%diag(zz$d[1:3])

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

# curves <- slingshot(naive, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
#                                 cluster_group_list = cluster_group_list)
# lineages <- .get_lineages(dat, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
#                           cluster_group_list = cluster_group_list)

#############

starting_cluster = cluster_group_list[[1]][1]
stopifnot(!is.list(cluster_group_list) || starting_cluster %in% cluster_group_list[[1]])
stopifnot(all(cluster_labels > 0), all(cluster_labels %% 1 == 0), length(unique(cluster_labels)) == max(cluster_labels))
if(all(!is.na(cluster_group_list))){
  tmp <- unlist(cluster_group_list)
  stopifnot(length(tmp) == length(unique(tmp)), length(tmp) == length(unique(cluster_labels)))
}

### construct the distance matrix
# dist_mat <- .compute_cluster_distances(dat, cluster_labels)

#############

k <- max(cluster_labels)
dist_mat <- matrix(0, k, k)

for(i in 1:(k-1)){
  for(j in (i+1):k){
    idx1 <- which(cluster_labels == i)
    idx2 <- which(cluster_labels == j)

    mean_vec1 <- colMeans(dat[idx1,])
    mean_vec2 <- colMeans(dat[idx2,])

    cov_mat1 <- stats::cov(dat[idx1,])
    cov_mat2 <- stats::cov(dat[idx2,])

    dist_mat[i,j] <- .covariance_distance(mean_vec1, cov_mat1, mean_vec2, cov_mat2)
    dist_mat[j,i] <- dist_mat[i,j]
  }
}

dist_mat
