# distribution: poisson
# natural parameter: m_{ij} = u_i^Tv_j
# relation to canonical parameters: m_{ij} = log(lambda)
# optimization problem: exp(m_{ij}) - a_{ij}*m_{ij}

.evaluate_objective.poisson <- function(dat, u_mat, v_mat){
  pred_mat <- u_mat %*% t(v_mat)
  idx <- which(!is.na(dat))
  stopifnot(all(pred_mat[idx] > 0))

  sum(exp(pred_mat)) - sum(dat[idx] * pred_mat[idx])
}

.evaluate_objective_single.poisson <- function(dat_vec, current_vec, other_mat){
  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))
  stopifnot(all(pred_vec[idx] > 0))

  sum(exp(pred_vec)) - sum(dat_vec[idx] * pred_vec[idx])
}

.gradient_vec.poisson <- function(dat_vec, current_vec, other_mat){
  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))

  tmp <- sapply(idx, function(j){
    other_mat[j,,drop=F]*(exp(pred_vec[j]) - dat_vec[j])
  })

  if(is.matrix(tmp)) rowSums(tmp) else sum(tmp)
}

