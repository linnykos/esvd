rm(list=ls())
load("../results/factorization_results.RData")

plot(obj$cell_mat[,1], obj$cell_mat[,2], asp = T,
     pch = 16, col = c(1:4)[rep(1:4, each = n_each)])

k <- 6
plot(res_our_list[[k]]$u_mat[,1], res_our_list[[k]]$u_mat[,2], asp = T,
     pch = 16, col = c(1:4)[rep(1:4, each = n_each)])

pred_mat <- 1/(res_our_NA_list[[k]]$u_mat %*% t(res_our_NA_list[[k]]$v_mat))
pred_mat <- t(sapply(1:nrow(pred_mat), function(x){
  pred_mat[x,] * extra_weight[x]
}))
plot(pred_mat[missing_idx], dat_impute[missing_idx], asp = T, pch = 16,
     col = rgb(0,0,0,0.2))
lines(c(0,1e5), c(0,1e5), col = "red", lwd = 2)



plot(res_ica[,1], res_ica[,2], asp = T,
     pch = 16, col = c(1:4)[rep(1:4, each = n_each)])

plot(res_svd[,1], res_svd[,2], asp = T,
     pch = 16, col = c(1:4)[rep(1:4, each = n_each)])
