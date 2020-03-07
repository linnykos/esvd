rm(list=ls())
load("../results/step5_clustering.RData")

set.seed(10)
# our_curves <- eSVD::slingshot(res_our$u_mat[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
#                               cluster_group_list = cluster_group_list,
#                               verbose = T, upscale_vec = upscale_vec,
#                               reduction_percentage = 0.25)
# our_curves$lineages

##########

dat <- res_our$u_mat[,1:p]
starting_cluster <- cluster_group_list[[1]][1]
# lineages <- .get_lineages(dat, cluster_labels, starting_cluster = starting_cluster,
#                           cluster_group_list = cluster_group_list,
#                           use_initialization = F)
# dist_mat <- .compute_cluster_distances(dat, cluster_labels)
# dist_mat = dist_mat^2
# dist_mat[c(2,3,11), c(2,3,11)]
#
# lineages <- .construct_lineage_from_hierarchy(dist_mat, cluster_group_list,
#                                               starting_cluster)
# lineages

.covariance_distance2 <- function(mean_vec1, cov_mat1, mean_vec2, cov_mat2, n1, n2, tol = 0.1){
  mat <- cov_mat1/n1 + cov_mat2/n2

  eigen_res <- eigen(mat)
  eigen_res$values[eigen_res$values < tol] <- tol
  mat <- eigen_res$vectors %*% diag(1/eigen_res$values) %*% t(eigen_res$vectors)

  as.numeric(t(mean_vec1 - mean_vec2) %*% mat %*% (mean_vec1 - mean_vec2))
}

.compute_cluster_distances2 <- function(dat, cluster_labels){
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

      dist_mat[i,j] <- .covariance_distance2(mean_vec1, cov_mat1, mean_vec2, cov_mat2,
                                             length(idx1), length(idx2))
      dist_mat[j,i] <- dist_mat[i,j]
    }
  }

  dist_mat
}

set.seed(10)
dist_mat2 <- .compute_cluster_distances2(dat, cluster_labels)
dist_mat2

lineages <- .construct_lineage_from_hierarchy(dist_mat2, cluster_group_list, starting_cluster)
lineages

