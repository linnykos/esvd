rm(list=ls())
load("../tmp.RData")
i = 4
tmp <- res_list[[i]]$u_mat; svd_res <- svd(tmp); tmp <- svd_res$u[,1:res$k] %*% diag(svd_res$d[1:res$k])
# curves <- slingshot(tmp, cluster_labels = rep(1:4, each = n_seq[i]),
#                     starting_cluster = 1)

###########

plot(tmp[,1], tmp[,2],
     pch = 16, col = col_vec[rep(1:4, each = n_seq[i])], asp = T,
     xlab = "X[,1]", ylab = "X[,2]")


lineages <- .get_lineages(tmp, cluster_labels = rep(1:4, each = n_seq[i]), starting_cluster = 1)

dat = tmp
cluster_labels = rep(1:4, each = n_seq[i])
starting_cluster = 1
knn = NA

cluster_mat <- .construct_cluster_matrix(cluster_labels)
k <- ncol(cluster_mat)
stopifnot(is.matrix(dat), nrow(dat) == nrow(cluster_mat))

### get the connectivity matrix
centers <- .compute_cluster_center(dat, cluster_mat)
dat_augment <- rbind(centers, dat)

### construct the k-nearest neighbor graph
if(is.na(knn)){
  knn <- 1
  while(TRUE){
    knn_graph <- .construct_knn_graph(dat_augment, knn = knn)
    if(igraph::components(knn_graph)$no == 1) break()
    knn <- knn + 1
  }
} else {
  knn_graph <- .construct_knn_graph(dat_augment, knn = knn)
  stopifnot(igraph::components(knn_graph)$no == 1)
}

# spt_graph <- .construct_spt(knn_graph, k = k, starting_cluster = starting_cluster)
dist_mat <- igraph::distances(knn_graph, v = 1:k, to = 1:k, mode = "all")
dist_mat[is.infinite(dist_mat)] <- 0
dist_mat <- dist_mat^2
