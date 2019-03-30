rm(list=ls())
set.seed(10)
dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
dat[sample(1:prod(dim(dat)), 10)] <- NA
u_mat <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
v_mat <- abs(matrix(rnorm(8), nrow = 4, ncol = 2))
pred_mat <- u_mat %*% t(v_mat)
class(dat) <- c("gaussian", class(dat)[length(class(dat))])
scalar = 2

idx <- which(!is.na(dat))

res1 <- sum(-log(pred_mat[idx]) -
      pred_mat[idx]*dat[idx]*scalar^2 +
      pred_mat[idx]^2*dat[idx]^2*scalar^2)

pred_mat2 <- pred_mat
############

extra_weights = rep(1, nrow(dat))
pred_mat <- u_mat %*% t(v_mat)
idx <- which(!is.na(dat))
stopifnot(all(pred_mat > 0))
extra_mat <- t(sapply(1:nrow(dat), function(x){rep(extra_weights[x], ncol(dat))}))

res2 <- sum(-log(pred_mat[idx]) -
      pred_mat[idx]*dat[idx]*scalar^2/extra_mat[idx] +
      pred_mat[idx]^2*dat[idx]^2*scalar^2/(2*extra_mat[idx]^2))
