# distribution: poisson
# natural parameter: m_{ij} = u_i^Tv_j
# relation to canonical parameters: m_{ij} = log(lambda_{ij})
# optimization problem: exp(m_{ij}) - a_{ij}*mu_{ij}

.evaluate_objective.poisson <- function(dat, u_mat, v_mat, ...){
  stopifnot(ncol(u_mat) == ncol(v_mat), nrow(u_mat) == nrow(dat), nrow(v_mat) == ncol(dat))

  n <- nrow(dat); p <- ncol(dat)
  pred_mat <- u_mat %*% t(v_mat)
  idx <- which(!is.na(dat))
  stopifnot(all(pred_mat[idx] > 0))
  stopifnot(all(dat[idx] >= 0))

  1/(n*p) * sum(exp(pred_mat[idx]) - pred_mat[idx]*dat[idx])
}

.evaluate_objective_single.poisson <- function(dat_vec, current_vec, other_mat, n, p, ...){
  stopifnot(length(current_vec) == ncol(other_mat), nrow(other_mat) == length(dat_vec))

  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))
  stopifnot(all(pred_vec[idx] > 0))
  stopifnot(all(dat_vec[idx] >= 0))

  1/(n*p) * sum(exp(pred_vec[idx]) - pred_vec[idx]*dat_vec[idx])
}

.gradient_vec.poisson <- function(dat_vec, current_vec, other_mat, n, p,  ...){
  stopifnot(length(current_vec) == ncol(other_mat), nrow(other_mat) == length(dat_vec))

  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))
  stopifnot(all(dat_vec[idx] >= 0))
  stopifnot(all(pred_vec[idx] > 0))

  tmp <- sapply(idx, function(j){
    other_mat[j,,drop=F]*(exp(pred_vec[j])-dat_vec[j])
  })

  if(is.matrix(tmp)) 1/(n*p) * rowSums(tmp) else 1/(n*p) * sum(tmp)
}

.evaluate_objective_mat.poisson <- function(dat, pred_mat, ...){
  stopifnot(all(dim(dat) == dim(pred_mat)))

  n <- nrow(dat); p <- ncol(dat)
  idx <- which(!is.na(dat))
  stopifnot(all(pred_mat > 0))
  stopifnot(all(dat[idx] >= 0))

  1/(n*p) * sum(exp(pred_mat[idx]) - pred_mat[idx]*dat[idx])
}

.gradient_mat.poisson <- function(dat, pred_mat, ...){
  stopifnot(all(dim(dat) == dim(pred_mat)))

  n <- nrow(dat); p <- ncol(dat)
  stopifnot(all(!is.na(dat)))
  stopifnot(all(pred_mat > 0))
  stopifnot(all(dat >= 0))

  (exp(pred_mat) - dat)/(n*p)
}



