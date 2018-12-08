rm(list=ls())
load("../results/factorization_results.RData")

cluster_labels <- c(1:4)[rep(1:4, each = paramMat[1,"n"])]
.b_estimate <- function(mat, cluster_labels){
  uniq_val <- sort(unique(cluster_labels))
  mean(sapply(uniq_val, function(x){
    idx <- which(cluster_labels == x)
    stats::princomp(mat[idx,])$sdev[1]
  }))/4
}

.compare_two_lineages <- function(target_lineage, source_lineage, cluster_labels){
  len1 <- length(target_lineage$curves)
  len2 <- length(source_lineage$curves)

  sapply(1:len1, function(x){
    target_lambda <- target_lineage$curves[[x]]$lambda_long
    idx <- which(target_lambda != 0)
    target_lambda <- target_lambda[idx]

    min(sapply(1:len2, function(y){
      source_lambda <- source_lineage$curves[[y]]$lambda_long
      source_lambda <- source_lambda[idx]
      abs(stats::cor(target_lambda, source_lambda))
    }))
  })
}

kendall_list <- vector("list", length(res[[1]]))
for(x in 1:length(kendall_list)){
  print(x)
  b_est <- .b_estimate(res[[1]][[x]]$cell_mat[,1:2], cluster_labels)
  true_lineage <- singlecell::slingshot(res[[1]][[x]]$cell_mat[,1:2], cluster_labels, 1, knn = NA,
                                        b = b_est, remove_outlier = F)
#
  plot(res[[1]][[x]]$cell_mat[,1], res[[1]][[x]]$cell_mat[,2], asp = T, pch = 16, col = rgb(0,0,0,0.1))
  for(i in 1:length(true_lineage$curves)){
    ord <- true_lineage$curves[[i]]$ord
    lines(true_lineage$curves[[i]]$s[ord, 1], true_lineage$curves[[i]]$s[ord, 2], lwd = 3,
          col = "black")
  }

  b_est <- .b_estimate(res[[1]][[x]]$res_our$u_mat[,1:2], cluster_labels)
  our_lineage <- singlecell::slingshot(res[[1]][[x]]$res_our$u_mat[,1:2], cluster_labels, 1, knn = NA,
                                       b = b_est)

  plot(res[[1]][[x]]$res_our$u_mat[,1], res[[1]][[x]]$res_our$u_mat[,2], asp = T, pch = 16, col = cluster_labels)
  for(i in 1:length(our_lineage$curves)){
    ord <- our_lineage$curves[[i]]$ord
    lines(our_lineage$curves[[i]]$s[ord, 1], our_lineage$curves[[i]]$s[ord, 2], lwd = 3,
          col = "black")
  }

  b_est <- .b_estimate(res[[1]][[x]]$res_svd[,1:2], cluster_labels)
  svd_lineage <- singlecell::slingshot(res[[1]][[x]]$res_svd[,1:2], cluster_labels, 1, knn = NA,
                                       b = b_est, percentage = 0.3)

  b_est <- .b_estimate(res[[1]][[x]]$res_ica[,1:2], cluster_labels)
  ica_lineage <- singlecell::slingshot(res[[1]][[x]]$res_ica[,1:2], cluster_labels, 1, knn = NA,
                                       b = b_est, percentage = 0.3)

  res_mat <- matrix(c(.compare_two_lineages(true_lineage, our_lineage),
                      .compare_two_lineages(true_lineage, svd_lineage),
                      .compare_two_lineages(true_lineage, ica_lineage)), ncol = 2, byrow = T)
  rownames(res_mat) <- c("Our", "SVD", "ICA")
  kendall_list[[x]] <- res_mat
}

our_vec <- unlist(lapply(kendall_list, function(x){x["Our",]}))
svd_vec <- unlist(lapply(kendall_list, function(x){x["SVD",]}))
ica_vec <- unlist(lapply(kendall_list, function(x){x["ICA",]}))

hist(our_vec, breaks = 25, col = "gray")
hist(svd_vec, breaks = 25, col = "gray")
hist(ica_vec, breaks = 25, col = "gray")
