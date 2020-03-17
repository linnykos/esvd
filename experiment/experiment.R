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

res <- tuning_scalar(dat, family = "neg_binom", search_min = 1, search_max = 2*max(dat),
                     max_val = 100, k = 1)


