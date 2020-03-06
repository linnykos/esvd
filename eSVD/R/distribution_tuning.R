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
tuning <- function(dat, u_mat, v_mat, family_to, family_from){
  if(family_to == "neg_binom" & family_from == "neg_binom"){
    .tuning_neg_binom(dat, u_mat, v_mat)
  } else if(family_to == "neg_binom" & family_from == "poisson"){
    .tuning_neg_binom_from_poisson(dat, u_mat, v_mat)
  } else if(family_to == "curved_gaussian" & family_from == "curved_gaussian"){
    .tuning_curved_gaussian(dat, u_mat, v_mat)
  } else {
    .tuning_curved_gaussian_from_exponential(dat, u_mat, v_mat)
  }
}

.tuning_neg_binom <- function(dat, u_mat, v_mat){
  p_mat <- exp(u_mat %*% t(v_mat))

  r_seq <- 1:100
  proposed_val <- t(sapply(r_seq, function(x){
    pred_mat <- x*p_mat/(1-p_mat)
    c(sum((dat - pred_mat)^2)/prod(dim(dat)),
      sum((pred_mat + pred_mat^2/x))/prod(dim(dat)))
  }))

  r_seq[which.min(abs(proposed_val[,1] - proposed_val[,2]))]
}

.tuning_neg_binom_from_poisson <- function(dat, u_mat, v_mat){
  pred_mat <- exp(u_mat %*% t(v_mat))

  target_val <- sum((dat - pred_mat)^2)/prod(dim(dat))
  r_seq <- 1:100
  proposed_val <- t(sapply(r_seq, function(x){
    sum((pred_mat + pred_mat^2/x))/prod(dim(dat))
  }))

  r_seq[which.min(abs(target_val - proposed_val))]
}

.tuning_curved_gaussian <- function(dat, u_mat, v_mat){
  pred_mat <- 1/(u_mat %*% t(v_mat))

  target_val <- sum((dat - pred_mat)^2)/prod(dim(dat))
  r_seq <- seq(1, 10, length.out = 101)
  proposed_val <- sapply(r_seq, function(x){sum((pred_mat/x)^2)/prod(dim(dat))})

  r_seq[which.min(abs(target_val - proposed_val))]
}

.tuning_curved_gaussian_from_exponential <- function(dat, u_mat, v_mat){
  pred_mat <- 1/(-u_mat %*% t(v_mat))

  target_val <- sum((dat - pred_mat)^2)/prod(dim(dat))
  r_seq <- seq(1, 10, length.out = 101)
  proposed_val <- sapply(r_seq, function(x){sum((pred_mat/x)^2)/prod(dim(dat))})

  r_seq[which.min(abs(target_val - proposed_val))]
}
