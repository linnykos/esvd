# distribution: one-parameter Gaussian where var = mean/2
# natural parameter: m_{ij} = u_i^Tv_j
# relation to canonical parameters: -2/mu^2 and m_{ij} = 4/mu
# optimization problem: -log(m_{ij}) - a_{ij}^2*(-m_{ij}^2/8) - a_{ij}*m_{ij}

.evaluate_objective.gaussian <- function(dat, u_mat, v_mat,
                                         extra_weights = rep(1, nrow(dat)), scalar = 2, ...){
  pred_mat <- u_mat %*% t(v_mat)
  idx <- which(!is.na(dat))
  stopifnot(all(pred_mat > 0))
  extra_mat <- t(sapply(1:nrow(dat), function(x){rep(extra_weights[x], ncol(dat))}))

  sum(-log(pred_mat[idx]) -
        pred_mat[idx]*dat[idx]*scalar^2/extra_mat[idx] +
        pred_mat[idx]^2*dat[idx]^2*scalar^2/(2*extra_mat[idx]^2))
}

.evaluate_objective_single.gaussian <- function(dat_vec, current_vec, other_mat,
                                                extra_weights = rep(1, nrow(other_mat)), scalar = 2, ...){
  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))
  stopifnot(all(pred_vec[idx] > 0))

  sum(-log(pred_vec[idx]) -
        pred_vec[idx]*dat_vec[idx]*scalar^2/extra_weights[idx] +
        pred_vec[idx]^2*dat_vec[idx]^2*scalar^2/(2*extra_weights[idx]^2))
}

.gradient_vec.gaussian <- function(dat_vec, current_vec, other_mat,
                                   extra_weights = rep(1, nrow(other_mat)), scalar = 2, ...){
  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))

  tmp <- sapply(idx, function(j){
    other_mat[j,,drop=F]*(pred_vec[j]*scalar^2*dat_vec[j]^2/extra_weights[j]^2 -
                            dat_vec[j]*scalar^2/extra_weights[j] - 1/pred_vec[j])
  })

  if(is.matrix(tmp)) rowSums(tmp) else sum(tmp)
}

.evaluate_objective_mat.gaussian <- function(dat, theta_mat, scalar = 2, ...){
  idx <- which(!is.na(dat))
  sum(-log(theta_mat[idx]) - (-scalar^2*dat[idx]*theta_mat[idx] + scalar^2*dat[idx]^2*theta_mat[idx]^2))
}

.gradient_mat.gaussian <- function(dat, theta_mat, scalar = 2, ...){
  idx <- which(!is.na(dat))
  sum(-1/(theta_mat[idx]) + scalar^2*dat[idx] - scalar^2*dat[idx]^2*theta_mat[idx])
}


