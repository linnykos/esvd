rm(list=ls())
load("../results/factorization_results.RData")

k <- 6
plot(res_our_list[[6]]$u_mat[,1], res_our_list[[6]]$u_mat[,2], asp = T)
