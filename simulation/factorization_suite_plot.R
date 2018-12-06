rm(list=ls())
load("../results/factorization_results.RData")

zz <- res[[1]][[1]]
plot(zz$res_our$u_mat[,1], zz$res_our$u_mat[,3], asp = T,
     pch = 16, col = c(1:4)[rep(1:4, each = paramMat[1,"n"])])

plot(zz$res_svd[,1], zz$res_svd[,2], asp = T,
     pch = 16, col = c(1:4)[rep(1:4, each = paramMat[1,"n"])])

plot(zz$res_ica[,1], zz$res_ica[,2], asp = T,
     pch = 16, col = c(1:4)[rep(1:4, each = paramMat[1,"n"])])

########

plot(obj$cell_mat[,1], obj$cell_mat[,2], asp = T,
     pch = 16, col = c(1:4)[rep(1:4, each = n_each)])

k <- 6


k <- 6
pred_mat <- 1/(res_our_NA_list[[k]]$u_mat %*% t(res_our_NA_list[[k]]$v_mat))
pred_mat <- t(sapply(1:nrow(pred_mat), function(x){
  pred_mat[x,] * extra_weight[x]
}))
plot(pred_mat[missing_idx], dat_impute[missing_idx], asp = T, pch = 16,
     col = rgb(0,0,0,0.2))
lines(c(0,1e5), c(0,1e5), col = "red", lwd = 2)

###########

.l2norm <- function(x){sqrt(sum(x^2))}
quality_vec <- sapply(res_our_NA_list, function(x){
  pred_mat <- 1/(x$u_mat %*% t(x$v_mat))
  pred_mat <- t(sapply(1:nrow(pred_mat), function(x){
    pred_mat[x,] * extra_weight[x]
  }))

  mat <- cbind(dat_impute[missing_idx], pred_mat[missing_idx])
  mat <- mat[which(mat[,1] >= 5), ]
  pca_res <- stats::princomp(mat)
  diag_vec <- c(1,1); diag_vec <- diag_vec/.l2norm(diag_vec)

  acos(diag_vec %*% pca_res$loadings[,1])
})



###########


plot(res_ica[,1], res_ica[,2], asp = T,
     pch = 16, col = c(1:4)[rep(1:4, each = n_each)])

plot(res_svd[,1], res_svd[,2], asp = T,
     pch = 16, col = c(1:4)[rep(1:4, each = n_each)])
