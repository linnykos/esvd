rm(list=ls())
x=1
set.seed(x*10)
dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))
class(dat) <- c("gaussian", class(dat))
bool <- sample(c(T, F), 1)

res <- .initialization(dat, family = "gaussian", max_val = 100)
u_mat <- res$u_mat
v_mat <- res$v_mat

# res <- .optimize_mat(dat, v_mat, u_mat, bool, max_val = 100)

left = bool
