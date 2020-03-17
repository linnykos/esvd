rm(list=ls())
set.seed(10)
n <- 150
u_mat <- abs(matrix(rnorm(n), nrow = n, ncol = 1))
v_mat <- -abs(matrix(rnorm(n), nrow = n, ncol = 1))
pred_mat <- u_mat %*% t(v_mat)
dat <- pred_mat

for(i in 1:n){
  for(j in 1:n){
    dat[i,j] <- stats::rnbinom(1, size = 50, prob = 1-exp(pred_mat[i,j]))
  }
}

res <- tuning_scalar(dat, family = "neg_binom", max_val = 100, k = 1)
param <- res[length(res)]

init <- initialization(dat, family = "neg_binom", scalar = param, max_val = 100, k = 1)
fit <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                  family =  "neg_binom", scalar = param,  max_val = 100, k = 1)

nat_mat <- fit$u_mat %*% t(fit$v_mat)
mean_mat <- compute_mean(nat_mat, family = "neg_binom", scalar = param)

plot(mean_mat, dat, asp = T)
