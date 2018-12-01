set.seed(10)
dat <- abs(matrix(rexp(20), nrow = 5, ncol = 4))

max_val <- 5
init <- .initialization(dat, family = "gaussian", scalar = 100, max_val = max_val)

stopifnot(all(init$u_mat %*% t(init$v_mat) <= max_val+1e-6))
quantile(init$u_mat %*% t(init$v_mat))

res <- .fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                          family = "gaussian",
                          max_val = max_val, scalar = 100)

stopifnot(all(res$u_mat %*% t(res$v_mat) <= max_val+1e-6))
