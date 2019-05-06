# distribution: exponential
# natural parameter: m_{ij} = u_i^Tv_j
# relation to canonical parameters: m_{ij} = -1/lambda
# optimization problem: -log(-m_{ij}) - a_{ij}*m_{ij}

.evaluate_objective.exponential <- function(dat, u_mat, v_mat, ...){
  n <- nrow(dat); p <- ncol(dat)
  pred_mat <- u_mat %*% t(v_mat)
  idx <- which(!is.na(dat))
  stopifnot(all(pred_mat[idx] < 0))

  1/(n*p) * sum(-log(-pred_mat[idx]) - pred_mat[idx]*dat[idx])
}

.evaluate_objective_single.exponential <- function(dat_vec, current_vec, other_mat, n, p, ...){
  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))
  stopifnot(all(pred_vec[idx] < 0))

  1/(n*p) * sum(-log(-pred_vec[idx]) - pred_vec[idx]*dat_vec[idx])
}

.gradient_vec.exponential <- function(dat_vec, current_vec, other_mat, n, p,  ...){
  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))

  tmp <- sapply(idx, function(j){
    other_mat[j,,drop=F]*(-1/pred_vec[j]-dat_vec[j])
  })

  if(is.matrix(tmp)) 1/(n*p) * rowSums(tmp) else 1/(n*p) * sum(tmp)
}


.evaluate_objective_mat.exponential <- function(dat, pred_mat, ...){
  n <- nrow(dat); p <- ncol(dat)
  stopifnot(all(pred_mat < 0))
  idx <- which(!is.na(dat))

  1/(n*p) * sum(-log(-pred_mat[idx]) - pred_mat[idx]*dat[idx])
}

.gradient_mat.exponential <- function(dat, pred_mat, ...){
  n <- nrow(dat); p <- ncol(dat)
  stopifnot(all(!is.na(dat)))

  (-1/(pred_mat) - dat)/(n*p)
}



