rm(list=ls())
load("../experiment/experiment.RData")

i <- 4
cluster_labels <- rep(1:4, each = vec["n_each"])
plot(fit_list[[i]]$u_mat[,1], fit_list[[i]]$u_mat[,2], pch = 16, col = cluster_labels, asp = T)

plot(fit$u_mat[,1], fit$u_mat[,2], pch = 16, col = cluster_labels, asp = T)
