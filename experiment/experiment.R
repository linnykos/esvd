rm(list=ls())
load("../results/factorization_results_others_tmp_gen4.RData")

cluster_labels <- rep(1:4, each = 50)
i <- 4
plot(res[[2]][[i]]$fit$u_mat[,1], res[[2]][[i]]$fit$u_mat[,2], asp = T, col = cluster_labels, pch = 16)

i <- 3
plot(res[[1]][[i]]$fit[,1], res[[1]][[i]]$fit[,2], asp = T, col = cluster_labels, pch = 16)

