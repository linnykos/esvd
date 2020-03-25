rm(list=ls())
load("../results/factorization_results_others_esvd.RData")

cluster_labels <- rep(1:4, each = 50)
i <- 1
par(mfrow = c(1,2))
plot(res[[1]][[i]]$fit$fit$u_mat[,1], res[[1]][[i]]$fit$fit$u_mat[,2], asp = T,
     col = cluster_labels, pch = 16, main = "Estimated")
plot(res[[1]][[i]]$truth[,1], res[[1]][[i]]$truth[,2], asp = T,
     col = cluster_labels, pch = 16, main = "Truth")


quality_vec <- sapply(1:length(res[[1]]), function(i){
  dist_mat_truth <- as.matrix(stats::dist(res[[1]][[i]]$truth))
  dist_mat_est <- as.matrix(stats::dist(res[[1]][[i]]$fit$fit$u_mat[,1:2]))

  mean(sapply(1:nrow(dist_mat_est), function(i){
    cor(dist_mat_truth[i,], dist_mat_est[i,], method = "kendall")
  }))
})


