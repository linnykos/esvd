rm(list=ls())
load("../results/step4_factorization.RData")

var_list <- ls()
var_list <- var_list[!var_list %in% c("dat_impute", "fitting_distr", "max_val")]
rm(list=var_list)

library(eSVD)
scalar <- 2
k <- 5

family = fitting_distr
max_val = max_val
max_iter = 10
tol = 1e-3
verbose = T
dat = dat_impute

session_info <- sessionInfo()

#############################

set.seed(10)
stopifnot(is.matrix(dat), k <= min(dim(dat)), k %% 1 == 0, k > 0)

direction <- eSVD:::.dictate_direction(family)
if(!is.na(max_val)){
  stopifnot(max_val > 0)

  # flip max_val so it's the correct sign downstream
  if(!is.na(direction) && direction == "<=") max_val <- -max_val
}

# initialize
dat2 <- eSVD:::.matrix_completion(dat, k = k)
if(length(class(dat2)) == 1) class(dat2) <- c(family, class(dat2)[length(class(dat2))])

# projected gradient descent
nat_mat <- eSVD:::.projected_gradient_descent(dat2, k = k, max_val = max_val,
                                       direction = direction,
                                       max_iter = max_iter,
                                       tol = tol, scalar = scalar)

res <- eSVD:::.svd_projection(nat_mat, k = k, factors = T)
u_mat <- res$u_mat; v_mat <- res$v_mat

if(!is.na(direction)){
  if(direction == "<=") {
    stopifnot(all(nat_mat[which(!is.na(dat))] < 0))
  } else {
    stopifnot(all(nat_mat[which(!is.na(dat))] > 0))
  }
}

tmp <- eSVD:::.fix_rank_defficiency_initialization(u_mat, v_mat, direction)
u_mat2 <- tmp$u_mat; v_mat2 <- tmp$v_mat

tmp <- eSVD:::.reparameterize(u_mat2, v_mat2)
u_mat3 <- tmp$u_mat; v_mat3 <- tmp$v_mat

rm(list = "tmp")

set.seed(10)
init <- eSVD::initialization(dat_impute, family = fitting_distr, k = k, max_val = max_val,
                             scalar = scalar)

save.image("../experiment/debugging_0037.RData")




