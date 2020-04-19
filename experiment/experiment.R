set.seed(10)
n <- 20
p <- 10
size <- 500
k <- 2
u_mat <- abs(matrix(rnorm(n), nrow = n, ncol = k))
v_mat <- -abs(matrix(rnorm(p), nrow = p, ncol = k))
pred_mat <- u_mat %*% t(v_mat)
dat <- pred_mat
class(dat) <- c("neg_binom", class(dat)[length(class(dat))])

for(i in 1:n){
  for(j in 1:p){
    dat[i,j] <- stats::rnbinom(1, size = size, prob = 1-exp(pred_mat[i,j]))
  }
}

missing_idx <- construct_missing_values(n = nrow(dat), p = ncol(dat), num_val = 1)
dat_NA <- dat
dat_NA[missing_idx] <- NA

scalar_vec <- c(10, 500, 10000)

fit_list <- lapply(scalar_vec, function(scalar){
  init <- initialization(dat_NA, family = "neg_binom", max_val = 100, scalar = scalar, k = k)
  fit <- fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                           max_iter = 10, max_val = 100, k = 1,
                           family = "neg_binom", scalar = scalar)
})

for(i in 1:3){
  plot_prediction_against_observed(dat, nat_mat_list = list(fit_list[[i]]$u_mat %*% t(fit_list[[i]]$v_mat)),
                                   family = "neg_binom", missing_idx_list = list(missing_idx),
                                   scalar = scalar_vec[i], plot = T, main = i)
}


#####

family = "neg_binom"
missing_idx_list = list(missing_idx)
width = 0.8

stopifnot(length(nat_mat_list_list) == length(scalar_vec))
stopifnot(length(unique(sapply(nat_mat_list_list, length))) == 1)
stopifnot(length(nat_mat_list_list[[1]]) == length(missing_idx_list))

nat_mat_list_list <- lapply(1:length(scalar_vec), function(i){
  list(fit_list[[i]]$u_mat %*% t(fit_list[[i]]$v_mat))
})

res_list <- lapply(1:length(nat_mat_list_list), function(i){
  plot_prediction_against_observed(dat, nat_mat_list_list[[i]], family = family,
                                   missing_idx_list = missing_idx_list,
                                   width = width, scalar = scalar_vec[i], plot = F)
})

