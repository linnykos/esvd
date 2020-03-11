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

missing_vec <- construct_missing_values(n = nrow(dat), p = ncol(dat), num_val = 2)
dat_NA <- dat
dat_NA[missing_vec] <- NA

#######

# fit initial fit
family_initial <- ifelse(family == "neg_binom", "poisson", "exponential")
fit <- .tuning_fit(dat_NA, family = family_initial, scalar = NA, max_val = 100, k = 1)

#######

# init <- initialization(dat_NA, family = family_initial, max_val = 100, k = 1)

#####

direction <- .dictate_direction(family)
if(!is.na(max_val)){
  if(!is.na(direction) && direction == "<=") max_val <- -max_val
}

# initialize
# dat <- .matrix_completion(dat_NA, k = k)

#############

