# distribution: one-parameter Gaussian where var = mean/scalar
# natural parameter: m_{ij} = u_i^Tv_j
# relation to canonical parameters: m_{ij} = 1/mu_{ij}
# optimization problem: -log(m_{ij}) - scalar^2*a_{ij}^2*(-m_{ij}^2)/n - scalar^2*a_{ij}*m_{ij}

.evaluate_objective.curved_gaussian <- function(dat, u_mat, v_mat, scalar = 2, ...){
  stopifnot(ncol(u_mat) == ncol(v_mat), nrow(u_mat) == nrow(dat), nrow(v_mat) == ncol(dat))

  n <- nrow(dat); p <- ncol(dat)
  pred_mat <- u_mat %*% t(v_mat)
  idx <- which(!is.na(dat))
  stopifnot(all(pred_mat > 0))

  1/(n*p) * sum(-log(pred_mat[idx]) -
        pred_mat[idx]*dat[idx]*scalar^2 +
        pred_mat[idx]^2*dat[idx]^2*scalar^2/2)
}

.evaluate_objective_single.curved_gaussian <- function(dat_vec, current_vec, other_mat, n, p,
                                                scalar = 2, ...){
  stopifnot(length(current_vec) == ncol(other_mat), nrow(other_mat) == length(dat_vec))

  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))
  stopifnot(all(pred_vec[idx] > 0))

  1/(n*p) * sum(-log(pred_vec[idx]) -
        pred_vec[idx]*dat_vec[idx]*scalar^2 +
        pred_vec[idx]^2*dat_vec[idx]^2*scalar^2/2)
}

.gradient_vec.curved_gaussian <- function(dat_vec, current_vec, other_mat, n, p, scalar = 2, ...){
  stopifnot(length(current_vec) == ncol(other_mat), nrow(other_mat) == length(dat_vec))

  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))
  stopifnot(all(pred_vec > 0))

  tmp <- sapply(idx, function(j){
    other_mat[j,,drop=F]*(- 1/pred_vec[j] - dat_vec[j]*scalar^2 +
                          pred_vec[j]*scalar^2*dat_vec[j]^2)
  })

  if(is.matrix(tmp)) 1/(n*p) * rowSums(tmp) else 1/(n*p) * sum(tmp)
}

.evaluate_objective_mat.curved_gaussian <- function(dat, pred_mat, scalar = 2, ...){
  stopifnot(all(dim(dat) == dim(pred_mat)))

  n <- nrow(dat); p <- ncol(dat)
  stopifnot(all(pred_mat > 0))
  idx <- which(!is.na(dat))

  1/(n*p) * sum(-log(pred_mat[idx]) -
        pred_mat[idx]*dat[idx]*scalar^2 +
        pred_mat[idx]^2*dat[idx]^2*scalar^2/2)
}

.gradient_mat.curved_gaussian <- function(dat, pred_mat, scalar = 2, ...){
  stopifnot(all(dim(dat) == dim(pred_mat)))

  n <- nrow(dat); p <- ncol(dat)
  stopifnot(all(!is.na(dat)))
  stopifnot(all(pred_mat > 0))

  (-1/(pred_mat) - scalar^2*dat + scalar^2*dat^2*pred_mat)/(n*p)
}


