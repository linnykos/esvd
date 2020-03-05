#' Tuning distributions
#'
#' @param dat data matrix
#' @param u_mat initial fit from \code{fit_factorization}
#' @param v_mat initial fit from \code{fit_factorization}
#' @param family character such as \code{"gaussian"} or \code{"exponential"}
#' @param r_min smallest value
#' @param r_max largest value
#'
#' @return value
#' @export
tuning <- function(dat, u_mat, v_mat, family, r_min = 1, r_max = 100){
  if(family == "neg_binom"){
    .tuning_neg_binom(dat, u_mat, v_mat, r_min, r_max)
  } else {
    .tuning_curved_gaussian(dat, u_mat, v_mat)
  }
}

.tuning_neg_binom <- function(dat, u_mat, v_mat, r_min = 1, r_max = 100){
  pred_mat <- exp(u_mat %*% t(v_mat))

  r_seq <- r_min:r_max
  vec <- sapply(r_seq, function(x){
    sum((dat - pred_mat)^2/(pred_mat + pred_mat^2/x))/2
  })

  max(r_seq[which.min(abs(vec - prod(dim(dat))))], 1)
}

.tuning_curved_gaussian <- function(dat, u_mat, v_mat){
  pred_mat <- -1/(u_mat %*% t(v_mat))

  max((prod(dim(dat))/(sum((dat - pred_mat)^2/(pred_mat)^2)))^(1/2), 1)
}
