.identification <- function(cov_x, cov_y, check = F){
  eigen_x <- eigen(cov_x)
  eigen_y <- eigen(cov_y)

  Vx <- eigen_x$vectors
  Vy <- eigen_y$vectors
  Dx <- diag(eigen_x$values);
  Dy <- diag(eigen_y$values)

  tmp <- sqrt(Dy) %*% t(Vy) %*% Vx %*% sqrt(Dx)
  svd_tmp <- svd(tmp)
  R <- svd_tmp$u %*% t(svd_tmp$v)

  # run a check
  if(check){
    Q <- t(R) %*% tmp
    stopifnot(sum(abs(Q - t(Q))) <= 1e-6)
  }

  Dx_inv <- Dx; diag(Dx_inv) <- 1/diag(Dx)
  sym_prod <- Vx %*% sqrt(Dx_inv) %*% t(R) %*% sqrt(Dy) %*% t(Vy)

  if(check){
    stopifnot(sum(abs(sym_prod - t(sym_prod))) <= 1e-6)
  }

  eigen_sym <- eigen(sym_prod)
  T_mat <- diag(sqrt(eigen_sym$values)) %*% t(eigen_sym$vectors)

  if(check){
    mat1 <- T_mat %*% cov_x %*% t(T_mat)

    T_mat_inv <- solve(T_mat)
    mat2 <- t(T_mat_inv) %*% cov_y %*% T_mat_inv
    stopifnot(sum(abs(mat1 - mat2)) <= 1e-6)
  }

  # adjust the transformation so it yields a diagonal matrix
  eig_res <- eigen(T_mat %*% cov_x %*% t(T_mat))

  t(eig_res$vectors) %*% T_mat
}

# cov_x = matrix(c(2,1,1,2),2,2); cov_y = matrix(c(5,-1,-1,5),2,2)
# T_mat = .identification(cov_x, cov_y)
# T_mat %*% cov_x %*% t(T_mat); solve(t(T_mat)) %*% cov_y %*% solve(T_mat)
