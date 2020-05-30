rm(list=ls())

set.seed(10)
n <- 20
p <- 10
k <- 2
xx <- matrix(rnorm(n*p), ncol = p, nrow = n)
svd_res <- svd(xx)

mat1 <- (n/p)^(1/4)*svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
mat2 <- (n/p)^(1/4)*xx %*% svd_res$v %*% diag(c(1/sqrt(svd_res$d[1:k]), rep(0, ncol(svd_res$v) - k)))

sum(abs(mat1 - mat2[,1:k]))
