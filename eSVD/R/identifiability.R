.identification <- function(cov_x, cov_y, check = F, tol = 1e-6){
  stopifnot(all(dim(cov_x) == dim(cov_y)), nrow(cov_x) == ncol(cov_x))
  if(nrow(cov_x) == 1){
    return(matrix((as.numeric(cov_y)/as.numeric(cov_x))^(1/4), 1, 1))
  }

  eigen_x <- eigen(cov_x)
  eigen_y <- eigen(cov_y)

  Vx <- eigen_x$vectors
  Vy <- eigen_y$vectors

  if(any(eigen_x$values <= tol) | any(eigen_y$values <= tol)) warning("Detecting rank defficiency in reparameterization step")

  Dx <- diag(eigen_x$values)
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


#' Function to reparameterize two matrices
#'
#' Designed to output matrices of the same dimension as \code{u_mat}
#' and \code{v_mat}, but linearly transformed so \code{u_mat %*% t(v_mat)}
#' is preserved but either \code{u_mat %*% t(u_mat)} is diagonal and equal to
#' \code{v_mat %*% t(v_mat)} (if \code{equal_covariance} is \code{FALSE})
#' or  \code{u_mat % *% t(u_mat)/nrow(u_mat)} is diagonal and equal to
#' \code{v_mat %*% t(v_mat)/nrow(v_mat)} (if \code{equal_covariance} is \code{TRUE})
#'
#' @param u_mat matrix of dimension \code{n} by \code{k}
#' @param v_mat matrix of dimension \code{p} by \code{k}
#' @param equal_covariance boolean
#'
#' @return list of two matrices
.reparameterize <- function(u_mat, v_mat, equal_covariance = F){
  stopifnot(ncol(u_mat) == ncol(v_mat))
  n <- nrow(u_mat); p <- nrow(v_mat)

  res <- .identification(t(u_mat) %*% u_mat, t(v_mat) %*% v_mat)

  if(equal_covariance){
    list(u_mat = (n/p)^(1/4)*u_mat %*% t(res), v_mat = (p/n)^(1/4)*v_mat %*% solve(res))
  } else {
    list(u_mat = u_mat %*% t(res), v_mat = v_mat %*% solve(res))
  }

}
