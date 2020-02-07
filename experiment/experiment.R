rm(list=ls())
x = 80
set.seed(10*x)
u_mat <- MASS::mvrnorm(5, rep(0, 5), diag(5))
v_mat <- MASS::mvrnorm(5, rep(1, 5), 2*diag(5))

res <- .reparameterize(u_mat, v_mat)

pred_mat <- u_mat %*% t(v_mat)
pred_mat2 <- res$u_mat %*% t(res$v_mat)

sum(abs(pred_mat - pred_mat2)) <= 1e-6

#####################

k <- ncol(u_mat)
pred_mat <- u_mat %*% t(v_mat)
svd_res <- svd::propack.svd(pred_mat, neig = k+1)

svd_res2 <- svd(pred_mat)

sum(abs(svd_res$d[1:k] - svd_res2$d[1:k]))
