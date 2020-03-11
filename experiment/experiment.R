set.seed(10)
u_mat <- abs(matrix(rnorm(60), nrow = 30, ncol = 2))
v_mat <- -abs(matrix(rnorm(60), nrow = 30, ncol = 2))
pred_mat <- u_mat %*% t(v_mat)
dat <- pred_mat

for(i in 1:30){
  for(j in 1:30){
    dat[i,j] <- stats::rnbinom(1, size = 50, prob = 1-exp(pred_mat[i,j]))
  }
}

init <- .tuning_fit(dat, family = "poisson", max_val = 100, max_iter = 10, k = 1)

search_min = 1
search_max = 2000
search_iter = 15
search_grid = 10

u_mat <- init$u_mat
v_mat <- init$v_mat
family <- "poisson"

#######

stopifnot(ncol(u_mat) == ncol(v_mat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat))
k <- ncol(u_mat); n <- nrow(dat); p <- ncol(dat)
df_val <- n*p - (n*k + p*k)
stopifnot(df_val > 0)

nat_mat <- u_mat %*% t(v_mat)
mean_mat_tmp <- compute_mean(nat_mat, family, scalar = 1)
recompute_mean <- family == "neg_binom"

lo_val <- search_min
hi_val <- search_max
min_val <- Inf
for(iter in 1:search_iter){
  scalar_seq <- seq(lo_val, hi_val, length.out = search_grid)
  obj_seq <- sapply(scalar_seq, function(scalar){
    .compute_tuning_objective(dat, family, nat_mat, mean_mat_tmp, scalar = scalar,
                              recompute_mean = recompute_mean)
  })
  plot(scalar_seq, obj_seq, main = iter)

  prev_min <- min_val
  min_val <- scalar_seq[which.min(abs(obj_seq - df_val))]

  #if(abs(prev_min - min_val) <= 1e-3) break()
  width <- abs(hi_val - lo_val)
  lo_val <- max(min_val - width/4, 1)
  hi_val <- min(min_val + width/4, search_max)
}
