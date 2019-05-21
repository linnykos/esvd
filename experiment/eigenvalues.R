rm(list=ls())
n_vec <- seq(10, 1000, length.out = 11)
vec <- sapply(1:length(n_vec), function(x){
  mat <- matrix(0, n_vec[x], n_vec[x])
  diag(mat) <- 1/n_vec[x]
  range(eigen(mat)$values)
})
