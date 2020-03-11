rm(list=ls())
set.seed(10)
scalar <- 50
prob <- 0.25
dat <- matrix(stats::rnbinom(200, size = scalar, prob = 1-prob), 10, 20)
# res <- tuning_scalar(dat, family = "neg_binom", max_val = 100, k = 1)

family = "neg_binom"
max_val = 100
k = 1
iter_max = 5
search_min = 1
search_max = 2000
search_iter = 10
search_grid = 10

#######

# fit initial fit
family_initial <- ifelse(family == "neg_binom", "poisson", "exponential")
fit <- .tuning_fit(dat, family = family_initial, scalar = NA, max_val = 100, k = 1)

###

# determine initial param
scalar_vec <- rep(NA, iter_max)
scalar_vec[1] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat, family = family_initial,
                                      search_min = search_min,
                                      search_max = search_max,
                                      search_iter = search_iter,
                                      search_grid = search_grid)

for(i in 2:iter_max){
  fit <- .tuning_fit(dat, family = family, scalar = scalar_vec[i-1], max_val = 100, k = 1)
  scalar_vec[i] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat, family = family,
                                        scalar = scalar_vec[i-1],
                                        search_min = search_min,
                                        search_max = search_max,
                                        search_iter = search_iter,
                                        search_grid = search_grid, max_val = 100, k = 1)
}
