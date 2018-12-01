# distribution: one-parameter Gaussian where var = mean/2
# natural parameter: m_{ij} = u_i^Tv_j
# relation to canonical parameters: -2/mu^2 and m_{ij} = 4/mu
# optimization problem: -log(m_{ij}) - a_{ij}^2*(-m_{ij}^2/8) - a_{ij}*m_{ij}

.evaluate_objective.gaussian <- function(dat, u_mat, v_mat,
                                         scalar = 2, ...){
  pred_mat <- u_mat %*% t(v_mat)
  idx <- which(!is.na(dat))
  if(any(pred_mat[idx] <= 0)){
    save(u_mat, v_mat, file = "../tmp.RData")
  }
  stopifnot(all(pred_mat[idx] > 0))

  sum(-log(pred_mat[idx]) - pred_mat[idx]*dat[idx] + pred_mat[idx]^2*dat[idx]^2/(2*scalar^2))
}

.evaluate_objective_single.gaussian <- function(dat_vec, current_vec, other_mat,
                                                scalar = 2, ...){
  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))
  stopifnot(all(pred_vec[idx] > 0))

  sum(-log(pred_vec[idx]) - pred_vec[idx]*dat_vec[idx] + pred_vec[idx]^2*dat_vec[idx]^2/(2*scalar^2))
}

.gradient_vec.gaussian <- function(dat_vec, current_vec, other_mat,
                                   scalar = 2, ...){
  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))

  tmp <- sapply(idx, function(j){
    other_mat[j,,drop=F]*(pred_vec[j]*dat_vec[j]^2/scalar^2 - dat_vec[j] - 1/pred_vec[j])
  })

  if(is.matrix(tmp)) rowSums(tmp) else sum(tmp)
}

