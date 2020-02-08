rm(list=ls())
set.seed(10)
n <- 100
theta_vec <- runif(n)
b_mat <- matrix(c(0.9, 0.2, 0.2, 0.9), 2, 2)
cluster_vec <- rep(1:2, each = n/2)

Matrix::rankMatrix(b_mat)

p_mat <- matrix(NA, n, n)
for(i in 1:n){
  for(j in 1:n){
    p_mat[i,j] <- theta_vec[i]*theta_vec[j]*b_mat[cluster_vec[i], cluster_vec[j]]
    p_mat[j,i] <- p_mat[i,j]
  }
}

Matrix::rankMatrix(p_mat)
