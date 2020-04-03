rm(list=ls())
load("../results/step5_clustering.RData")

.compute_cluster_distances_tmp <- function(dat, cluster_labels){
  k <- max(cluster_labels)
  dist_mat <- matrix(0, k, k)

  for(i in 1:(k-1)){
    for(j in (i+1):k){
      idx1 <- which(cluster_labels == i)
      idx2 <- which(cluster_labels == j)

      dat1 <- dat[idx1,]; dat2 <- dat[idx2,]
      if(nrow(dat1) < nrow(dat2)){
        dat1 <- dat1[sample(1:nrow(dat1), nrow(dat2), replace = T),]
      } else {
        dat2 <- dat2[sample(1:nrow(dat2), nrow(dat1), replace = T),]
      }

      dist_mat[i,j] <- transport::wasserstein(transport::pp(dat1),
                                              transport::pp(dat2), p = 1)

      # dist_mat[i,j] <- .covariance_distance_tmp(mean_vec1, cov_mat1, length(idx1),
      #                                       mean_vec2, cov_mat2, length(idx2))
      dist_mat[j,i] <- dist_mat[i,j]
    }
  }

  dist_mat
}

# .covariance_distance_tmp <- function(mean_vec1, cov_mat1, n1, mean_vec2, cov_mat2, n2, tol = 0.1){
#   mat <- cov_mat1/n1 + cov_mat2/n2
#
#   eigen_res <- eigen(mat)
#   eigen_res$values[eigen_res$values < tol] <- tol
#   mat <- eigen_res$vectors %*% diag(1/eigen_res$values) %*% t(eigen_res$vectors)
#
#   as.numeric(t(mean_vec1 - mean_vec2) %*% mat %*% (mean_vec1 - mean_vec2))
# }


#################################

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

upscale_vec <- rep(NA, length(unique(cluster_labels)))
size_vec <- sapply(cluster_group_list, function(x){length(which(cluster_labels %in% x))})
for(i in 1:length(cluster_group_list)){
  upscale_vec[cluster_group_list[[i]]] <- (max(size_vec)/size_vec[i])^(1/2)
}

p <- 5
dat <- esvd_embedding$u_mat[,1:p]
starting_cluster = cluster_group_list[[1]][1]
verbose = T
reduction_percentage = 0.25
use_initialization = F

########

stopifnot(!is.list(cluster_group_list) || starting_cluster %in% cluster_group_list[[1]])
stopifnot(all(cluster_labels > 0), all(cluster_labels %% 1 == 0), length(unique(cluster_labels)) == max(cluster_labels))
if(all(!is.na(cluster_group_list))){
  tmp <- unlist(cluster_group_list)
  stopifnot(length(tmp) == length(unique(tmp)), length(tmp) == length(unique(cluster_labels)))
}

### construct the distance matrix
dist_mat <- .compute_cluster_distances_tmp(dat, cluster_labels)
# dist_mat
# dist_mat <- dist_mat^2
zz <- dist_mat[c(3:9), c(3:9)]
colnames(zz) <- 3:9
rownames(zz) <- 3:9
zz

lineages <- .construct_lineage_from_hierarchy(dist_mat, cluster_group_list,
                                                starting_cluster)

lineages
