rm(list=ls())
load("../results/factorization_results_others.RData")

cluster_labels <- rep(1:4, each = 50)
i <- 1
par(mfrow = c(1,2))
plot(res[[1]][[i]]$fit$fit$u_mat[,1], res[[1]][[i]]$fit$fit$u_mat[,2], asp = T,
     col = cluster_labels, pch = 16, main = "Estimated")
plot(res[[1]][[i]]$truth[,1], res[[1]][[i]]$truth[,2], asp = T,
     col = cluster_labels, pch = 16, main = "Truth")


