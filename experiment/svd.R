rm(list=ls())

load("Zeisel_expr.RData")
#library(irlba)

dat <- as.matrix(dat)
#res <- irlba::irlba(dat)

res <- svd(dat)
save(res, "svd_tmp.RData")

res2 <- res
res2$d <- res$d[1:6]
res2$v <- res$v[,1:6]
res2$u <- res$u[,1:6]

dat_approx <- res2$u %*% diag(res2$d) %*% t(res2$v)
residual <- dat - dat_approx
