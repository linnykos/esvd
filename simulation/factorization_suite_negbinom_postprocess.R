rm(list=ls())
load("../results/factorization_results_negbinom_esvd.RData")

cluster_labels <- rep(1:4, each = 50)
j <- 1
i <- 32
par(mfrow = c(1,2))
plot(res[[j]][[i]]$fit$fit$u_mat[,1], res[[j]][[i]]$fit$fit$u_mat[,2], asp = T,
     col = cluster_labels, pch = 16, main = "Estimated")
plot(res[[j]][[i]]$truth[,1], res[[j]][[i]]$truth[,2], asp = T,
     col = cluster_labels, pch = 16, main = "Truth")

quality_vec <- sapply(1:length(res[[1]]), function(i){
  dist_mat_truth <- as.matrix(stats::dist(res[[1]][[i]]$truth))
  dist_mat_est <- as.matrix(stats::dist(res[[1]][[i]]$fit$fit$u_mat[,1:2]))

  mean(sapply(1:nrow(dist_mat_est), function(i){
    cor(dist_mat_truth[i,], dist_mat_est[i,], method = "kendall")
  }))
})
