.identification <- function(X, Y){
  stopifnot(ncol(X) == ncol(Y))

  cov_x <- t(X) %*% X
  cov_y <- t(Y) %*% Y

  eigen_x <- eigen(cov_x)
  eigen_y <- eigen(cov_y)

  Ux <- eigen_x$vectors
  Uy <- eigen_y$vectors
  Dx <- diag(eigen_x$values)
  Dy <- diag(eigen_y$values)

  inner_left <- sqrt(Dx) %*% t(Ux) %*% Uy %*% sqrt(Dy)
  svd_res <- svd(inner_left)
  Q <- svd_res$v %*% t(svd_res$u)

  sym_prod <- Ux %*% t(Ux) %*% Uy %*% sqrt(Dy) %*% Q %*% diag(1/sqrt(diag(Dx))) %*% t(Ux)

  eigen_sym <- eigen(sym_prod)
  A <- diag(sqrt(eigen_sym$values)) %*% t(eigen_sym$vectors)

  X <- X %*% t(A); Y <- Y %*% solve(A)

  list(X = X, Y = Y)
}
