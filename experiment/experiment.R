rm(list=ls())
set.seed(10)
scalar <- 50
prob <- 0.25
dat <- matrix(stats::rnbinom(200, size = scalar, prob = 1-prob), 10, 20)
u_mat <- matrix(1, nrow = nrow(dat), ncol = 1)
v_mat <- matrix(log(prob), nrow = ncol(dat), ncol = 1)
family = "neg_binom"

nat_mat <- u_mat %*% t(v_mat)
k <- ncol(u_mat); n <- nrow(dat); p <- ncol(dat)
df_val <- n*p - (n*k + p*k)

fn <- function(x){
  mean_mat <- compute_mean(nat_mat, family, scalar = x)
  var_mat <- .compute_variance(mean_mat, family, scalar = x)
  abs(sum((dat-mean_mat)^2/var_mat) - df_val)
}

gr <- function(x){
  mean_mat <- compute_mean(nat_mat, family, scalar = x)
  var_mat <- .compute_variance(mean_mat, family, scalar = x)
  .compute_gradient(dat, mean_mat, var_mat, df_val, family, scalar = x)
}

search_min <- 1
search_max <- 200

res <- stats::optim(par = (search_min+search_max)/2, fn = fn, gr = gr,
                    method = "L-BFGS-B",
                    lower = search_min, upper = search_max)
