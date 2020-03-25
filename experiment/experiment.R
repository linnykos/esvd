rm(list=ls())
x <- 46
print(x)
set.seed(10*x)
u_mat <- MASS::mvrnorm(10, rep(0, 5), diag(5))
v_mat <- MASS::mvrnorm(10, rep(1, 5), 2*diag(5))

res <- .reparameterize(u_mat, v_mat)

cov_x <- t(u_mat) %*% u_mat
cov_y <- t(v_mat) %*% v_mat

