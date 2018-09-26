set.seed(10)

X <- MASS::mvrnorm(n = 100, mu = rep(0,5), Sigma = diag(5))
Y <- MASS::mvrnorm(n = 100, mu = rep(1,5), Sigma = 2*diag(5))

Cx <- t(X)%*%X
Cy <- t(Y)%*%Y

eigen_x <- eigen(Cx)
eigen_y <- eigen(Cy)

tmp <- t(eigen_x$vectors)%*%eigen_y$vectors %*% diag(sqrt(eigen_y$values)) %*% diag(1/sqrt(eigen_x$values))
Matrix::rankMatrix(tmp)

eigen_tmp <- eigen(tmp)
sqrt_tmp <- eigen_tmp$vectors %*% sqrt(eigen) #uh-oh

tmp2 <- t(eigen_x$vectors)%*%eigen_y$vectors
eigen(tmp2)

### do eigenspaces commute?
zz <- t(eigen_x$vectors)%*%eigen_y$vectors
zz2 <- eigen_y$vectors%*%t(eigen_x$vectors)
