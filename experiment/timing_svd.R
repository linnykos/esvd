rm(list=ls())

d <- 200
k <- 2
mat1 <- matrix(rnorm(d*k), d, k)
mat <- mat1 %*% t(mat1)
microbenchmark::microbenchmark(svd(mat), svd::propack.svd(mat, neig = k))
