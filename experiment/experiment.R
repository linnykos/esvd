rm(list=ls())
set.seed(10)
n <- 100
u_mat <- abs(matrix(rnorm(2*n), nrow = n, ncol = 2))
v_mat <- -abs(matrix(rnorm(2*n), nrow = n, ncol = 2))
pred_mat <- u_mat %*% t(v_mat)
dat <- pred_mat

for(i in 1:n){
  for(j in 1:n){
    dat[i,j] <- stats::rnbinom(1, size = 50, prob = 1-exp(pred_mat[i,j]))
  }
}

family = "neg_binom"
max_val = 100
k = 1
iter_max = 5
search_min = 1
search_max = 2000
search_iter = 10
search_grid = 10

missing_vec <- construct_missing_values(n = nrow(dat), p = ncol(dat), num_val = 2)
dat_NA <- dat
dat_NA[missing_vec] <- NA

#######

# fit initial fit
family_initial <- ifelse(family == "neg_binom", "poisson", "exponential")
fit <- .tuning_fit(dat_NA, family = family_initial, scalar = NA, max_val = 100, k = 1)

###

# determine initial param
scalar_vec <- rep(NA, iter_max)
missing_idx <- which(is.na(dat))
scalar_vec[1] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat, family = family_initial,
                                      search_min = search_min,
                                      search_max = search_max,
                                      search_iter = search_iter,
                                      search_grid = search_grid,
                                      idx = missing_idx)

for(i in 2:iter_max){
  fit <- .tuning_fit(dat_NA, family = family, scalar = scalar_vec[i-1], max_val = 100, k = 1)
  scalar_vec[i] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat, family = family,
                                        scalar = scalar_vec[i-1],
                                        search_min = search_min,
                                        search_max = search_max,
                                        search_iter = search_iter,
                                        search_grid = search_grid, max_val = 100, k = 1,
                                        idx = missing_idx)
}

nat_mat <- fit$u_mat %*% t(fit$v_mat)

plot(pred_mat, nat_mat, asp = T, pch = 16)
lines(c(-1e5,1e5), c(-1e5,1e5), col = "red", lwd = 2, lty = 2)
