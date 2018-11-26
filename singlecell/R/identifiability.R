.identification <- function(cov_x, cov_y, check = T){
  eigen_x <- eigen(cov_x)
  eigen_y <- eigen(cov_y)

  Vx <- eigen_x$vectors
  Vy <- eigen_y$vectors
  Dx <- diag(eigen_x$values); Dx_inv <- Dx; diag(Dx_inv) <- 1/diag(Dx)
  Dy <- diag(eigen_y$values)

  tmp <- t(Vx)%*%Vy%*%sqrt(Dy)
  svd_tmp <- svd(tmp)
  R <- svd_tmp$u %*% diag(svd_tmp$d) %*% t(svd_tmp$u) %*% sqrt(Dx_inv)

  # run a check
  if(check){
    Q <- svd_tmp$u %*% t(svd_tmp$v)
    R2 <- t(Vx)%*%Vy%*%sqrt(Dy) %*% t(Q) %*% sqrt(Dx_inv)
    stopifnot(sum(abs(R2 - R)) <= 1e-6)
  }

  sym_prod <- Vx %*% R %*% t(Vx)

  if(check){
    stopifnot(sum(abs(sym_prod - t(sym_prod))) <= 1e-6)
  }

  eigen_sym <- eigen(sym_prod)
  T_mat <- diag(sqrt(eigen_sym$values)) %*% t(eigen_sym$vectors)

  if(check){
    mat1 <- t(T_mat)%*%T_mat%*%Vx%*%sqrt(Dx)%*%Q
    mat2 <- Vy %*% sqrt(Dy)
    stopifnot(sum(abs(mat1 - mat2)) <= 1e-6)
  }

  T_mat
}

# cov_x = matrix(c(2,1,1,2),2,2); cov_y = matrix(c(5,-1,-1,5),2,2)
# T_mat = .identification(cov_x, cov_y)
# T_mat %*% cov_x %*% t(T_mat); solve(t(T_mat)) %*% cov_y %*% solve(T_mat)
