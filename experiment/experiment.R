rm(list=ls())
p = 10

# generate the eigenvalues
set.seed(10)
mat <- matrix(rnorm(p^2), p, p)
mat <- mat %*% t(mat)
eig <- eigen(mat)
vec <- eig$vectors

target <- matrix(0, p, p)
for(i in 1:p){
  target <- target + 1/p * #(runif(1, 1, 2)) *
    vec[i,] %*% t(vec[i,])
}

eigen(target)$values

###################

mat <- diag(5)
mat[2,2] <- 0.3
eigen(mat)
