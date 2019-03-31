rm(list=ls())
set.seed(10)
dat <- abs(matrix(rexp(25, 1/2), nrow = 5, ncol = 5))
class(dat) <- c("exponential", class(dat)[length(class(dat))])
init <- initialization(dat, family = "exponential")

fit <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                         max_iter = 5, max_val = -100,
                         family = "exponential")


