rm(list=ls())
set.seed(10)
u_mat <- abs(matrix(rnorm(60)^2, nrow = 30, ncol = 2))
v_mat <- -abs(matrix(rnorm(60)^2, nrow = 30, ncol = 2))
pred_mat <- u_mat %*% t(v_mat)
dat <- pred_mat

for(i in 1:30){
  for(j in 1:30){
    dat[i,j] <- stats::rnbinom(1, size = 50, prob = 1-exp(pred_mat[i,j]))
  }
}

family = "neg_binom"
iter_max = 5
search_min = 1
search_max = 2*max(dat)
search_iter = 10
search_grid = 10
verbose = F

# fit initial fit
family_initial <- ifelse(family == "neg_binom", "poisson", "exponential")
fit <- .tuning_fit(dat, family = family_initial, scalar = NA, max_val = 100, k = 1)
if(verbose) print("Finished initial fit")

# determine initial param
scalar_vec <- rep(NA, iter_max)
scalar_vec[1] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat, family = family_initial,
                                      search_min = search_min,
                                      search_max = search_max, max_val = 100, k = 1)


# iterate between fitting and parameter estimation
for(i in 2:iter_max){
  fit <- .tuning_fit(dat, family = family, scalar = scalar_vec[i-1], max_val = 100, k = 1)
  if(verbose) print(paste0("Finished fit on iteration ", i))

  scalar_vec[i] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat, family = family,
                                        scalar = scalar_vec[i-1],
                                        search_min = search_min,
                                        search_max = search_max, max_val = 100, k = 1)

  if(verbose) print(paste0("Finished search on iteration ", i))
}

######

u_mat = fit$u_mat
v_mat = fit$v_mat
#family = "neg_binom"
family = family_initial

stopifnot(ncol(u_mat) == ncol(v_mat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat))
k <- ncol(u_mat); n <- nrow(dat); p <- ncol(dat)
df_val <- n*p - (n*k + p*k)

stopifnot(df_val > 0)

nat_mat <- u_mat %*% t(v_mat)

fn <- function(x){
  mean_mat <- compute_mean(nat_mat, family, scalar = x)
  var_mat <- .compute_variance(mean_mat, family, scalar = x)
  abs(sum((dat-mean_mat)^2/var_mat) - df_val)
}

fn(scalar_vec[1])
fn(50)
