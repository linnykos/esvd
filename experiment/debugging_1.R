rm(list=ls())
load("../results/factorization_results.RData")

x = 2
cluster_labels <- c(1:4)[rep(1:4, each = paramMat[1,"n"])]
b_est <- .b_estimate(res[[1]][[x]]$res_svd[,1:2], cluster_labels)
# svd_lineage <- singlecell::slingshot(res[[1]][[x]]$res_svd[,1:2], cluster_labels, 1, knn = NA,
#                                      b = b_est)

#######################

dat = res[[1]][[x]]$res_svd[,1:2]
starting_cluster = 1
knn = NA
remove_outlier = T
percentage = 0.05
shrink = 1
thresh = 0.001
max_iter = 15
b = b_est
upscale_vec = NA

lineages <- .get_lineages(dat, cluster_labels, starting_cluster = starting_cluster,
                          knn = knn, remove_outlier = remove_outlier,
                          percentage = percentage)


if(!any(is.na(upscale_vec))){
  idx_all <- unlist(lapply(1:max(cluster_labels, na.rm = T), function(x){
    idx <- which(cluster_labels == x)
    sample(idx, upscale_vec[x]*length(idx), replace = T)
  }))
  dat <- dat[idx_all,]
  cluster_labels <- cluster_labels[idx_all]
}

### setup
num_lineage <- length(lineages)
if(any(is.na(cluster_labels))) {
  idx <- which(is.na(cluster_labels))
  dat <- dat[-idx,]
  cluster_labels <- cluster_labels[-idx]
}
cluster_mat <- .construct_cluster_matrix(cluster_labels)
cluster_vec <- 1:ncol(cluster_mat)
centers <- .compute_cluster_center(dat, cluster_mat)

W <- .initialize_weight_matrix(cluster_mat, lineages)

### initial curves are piecewise linear paths through the tree
s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)
pcurve_list <- res$pcurve_list; D <- res$D

