rm(list=ls())
set.seed(10)

X <- MASS::mvrnorm(n = 100, mu = rep(0,5), Sigma = diag(5))
Y <- MASS::mvrnorm(n = 100, mu = rep(1,5), Sigma = 2*diag(5))

Cx <- t(X)%*%X
Cy <- t(Y)%*%Y

eigen_x <- eigen(Cx)
eigen_y <- eigen(Cy)

Ux <- eigen_x$vectors
Uy <- eigen_y$vectors
Dx <- diag(eigen_x$values)
Dy <- diag(eigen_y$values)

##########

inner_left <- sqrt(Dx) %*% t(Ux) %*% Uy %*% sqrt(Dy)
svd_res <- svd(inner_left)
Q <- svd_res$v %*% t(svd_res$u)

# check 1:
sum(abs(inner_left %*% Q - t(inner_left %*% Q))) <= 1e-6
# check 2:
tmp <- t(Ux) %*% Uy %*% sqrt(Dy) %*% Q %*% diag(1/sqrt(diag(Dx)))
sum(abs(tmp - t(tmp))) <= 1e-6

sym_prod <- Ux %*% t(Ux) %*% Uy %*% sqrt(Dy) %*% Q %*% diag(1/sqrt(diag(Dx))) %*% t(Ux)
#check 3:
sum(abs(sym_prod - t(sym_prod))) <= 1e-6

eigen_sym <- eigen(sym_prod)
A <- diag(sqrt(eigen_sym$values)) %*% t(eigen_sym$vectors)

# verify it works
val1 <- t(X[1,]) %*% Y[1,]
val2 <- t(A %*% X[1,]) %*% solve(t(A)) %*% Y[1,]

AX <- t(apply(X, 1, function(x){A%*%x}))
AinvtY <- t(apply(Y, 1, function(x){solve(t(A))%*%x}))

Cax <- t(AX) %*% AX
Cay <- t(AinvtY) %*% AinvtY
