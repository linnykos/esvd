set.seed(10)
n <- 100
u_mat <- abs(matrix(rnorm(n), nrow = n, ncol = 1))
v_mat <- -abs(matrix(rnorm(n), nrow = n, ncol = 1))
pred_mat <- u_mat %*% t(v_mat)
dat <- pred_mat

for(i in 1:n){
  for(j in 1:n){
    dat[i,j] <- stats::rnbinom(1, size = 50, prob = 1-exp(pred_mat[i,j]))
  }
}

mean(dat)
var(as.numeric(dat))
# res <- tuning_scalar(dat, family = "neg_binom", max_val = 100, k = 1)

family = "neg_binom"
iter_max = 5
search_min = 1
search_max = 2000
verbose = F

stopifnot(search_max > search_min)
missing_idx <- eSVD::construct_missing_values(n = nrow(dat), p = ncol(dat))

# fit initial fit
family_initial <- ifelse(family == "neg_binom", "poisson", "exponential")
dat_NA <- dat; dat_NA[missing_idx] <- NA
missing_val <- dat[missing_idx]
fit <- .tuning_fit(dat_NA, family = family_initial, scalar = NA,  max_val = 100, k = 1)
if(verbose) print("Finished initial fit")

# determine initial param
scalar_vec <- rep(NA, iter_max)
scalar_vec[1] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat,
                                      family = family_initial,
                                      idx = missing_idx,
                                      search_min = search_min,
                                      search_max = search_max,  max_val = 100, k = 1)


# iterate between fitting and parameter estimation
for(i in 2:iter_max){
  fit <- .tuning_fit(dat_NA, family = family, scalar = scalar_vec[i-1], max_val = 100, k = 1)
  if(verbose) print(paste0("Finished fit on iteration ", i))

  scalar_vec[i] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat, family = family,
                                        scalar = scalar_vec[i-1],
                                        idx = missing_idx,
                                        search_min = search_min,
                                        search_max = search_max, max_val = 100, k = 1)

  if(verbose) print(paste0("Finished search on iteration ", i))
}

################

u_mat <- fit$u_mat
v_mat <- fit$v_mat
idx <- missing_idx
# family <- family_initial
family <- "neg_binom"

stopifnot(ncol(u_mat) == ncol(v_mat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat))
k <- ncol(u_mat); n <- nrow(dat); p <- ncol(dat)
if(length(idx) == prod(dim(dat))){
  df_val <- n*p - (n*k + p*k)
} else {
  df_val <- length(idx)
}

stopifnot(df_val > 0)

nat_mat <- u_mat %*% t(v_mat)


fn <- function(x){
  mean_mat <- compute_mean(nat_mat, family, scalar = x)
  if(family %in% c("neg_binom")){
    abs(sum(dat[idx]/mean_mat[idx]) - df_val)
  } else {
    var_mat <- .compute_variance(mean_mat, family, scalar = x)
    abs(sum((dat[idx]-mean_mat[idx])^2/var_mat[idx]) - df_val)
  }
}

res <- stats::optimize(fn, interval = c(search_min, search_max))

fn(1)
fn(50)

mean_mat <- compute_mean(nat_mat, family, scalar = 1)
mean_mat[1:5,1:5]
dat[1:5,1:5]

##

x <- scalar_vec[length(scalar_vec)]
mean_mat <- compute_mean(nat_mat, family, scalar = x)
var_mat <-  .compute_variance(mean_mat, family, scalar = x)
abs(sum((dat[idx]-mean_mat[idx])^2/var_mat[idx]) - df_val)

plot(mean_mat[idx], dat[idx], asp = T)
