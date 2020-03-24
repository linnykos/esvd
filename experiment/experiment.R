rm(list=ls())
load("../results/factorization_results_others.RData")

cluster_labels <- rep(1:4, each = 50)
i <- 1
plot(res[[2]][[i]]$fit$fit$u_mat[,1], res[[2]][[i]]$fit$fit$u_mat[,2], asp = T, col = cluster_labels, pch = 16)

i <- 1
plot(res[[1]][[i]]$fit[,1], res[[1]][[i]]$fit[,2], asp = T, col = cluster_labels, pch = 16)

