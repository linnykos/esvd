# distribution: one-parameter Gaussian where var = mean/2
# natural parameter: m_{ij} = u_i^Tv_j
# relation to canonical parameters: -2/mu^2 and m_{ij} = 4/mu
# optimization problem: -log(m_{ij}) - a_{ij}^2*(-m_{ij}^2/8) - a_{ij}*m_{ij}

.evaluate_objective.gaussian <- function(dat, u_mat, v_mat, scalar = 2, ...){
  n <- nrow(dat); p <- ncol(dat)
  pred_mat <- u_mat %*% t(v_mat)
  idx <- which(!is.na(dat))
  stopifnot(all(pred_mat > 0))

  1/(n*p) * sum(-log(pred_mat[idx]) -
        pred_mat[idx]*dat[idx]*scalar^2 +
        pred_mat[idx]^2*dat[idx]^2*scalar^2/2)
}

.evaluate_objective_single.gaussian <- function(dat_vec, current_vec, other_mat, n, p,
                                                scalar = 2, ...){
  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))
  stopifnot(all(pred_vec[idx] > 0))

  1/(n*p) * sum(-log(pred_vec[idx]) -
        pred_vec[idx]*dat_vec[idx]*scalar^2 +
        pred_vec[idx]^2*dat_vec[idx]^2*scalar^2/2)
}

.gradient_vec.gaussian <- function(dat_vec, current_vec, other_mat, n, p, scalar = 2, ...){
  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))

  tmp <- sapply(idx, function(j){
    other_mat[j,,drop=F]*(pred_vec[j]*scalar^2*dat_vec[j]^2 -
                            dat_vec[j]*scalar^2 - 1/pred_vec[j])
  })

  if(is.matrix(tmp)) 1/(n*p) * rowSums(tmp) else 1/(n*p) * sum(tmp)
}

.evaluate_objective_mat.gaussian <- function(dat, pred_mat, scalar = 2, ...){
  n <- nrow(dat); p <- ncol(dat)
  stopifnot(all(pred_mat > 0))
  idx <- which(!is.na(dat))

  1/(n*p) * sum(-log(pred_mat[idx]) -
        pred_mat[idx]*dat[idx]*scalar^2 +
        pred_mat[idx]^2*dat[idx]^2*scalar^2/2)
}

.gradient_mat.gaussian <- function(dat, pred_mat, scalar = 2, ...){
  n <- nrow(dat); p <- ncol(dat)
  stopifnot(all(!is.na(dat)))

  (-1/(pred_mat) - scalar^2*dat + scalar^2*dat^2*pred_mat)/(n*p)
}


