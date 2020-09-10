#' Compute mean from natural parameter matrix
#'
#' @param nat_mat matrix
#' @param family character
#' @param ... additional parameters for distribution
#'
#' @return matrix
#' @export
compute_mean <- function(nat_mat, family, ...){
  if(family == "gaussian" || family == "poisson") {
    return(nat_mat)
  } else if(family == "neg_binom"){
    stopifnot(all(nat_mat < 0))
    .compute_mean_neg_binom(nat_mat, ...)
  } else if(family == "exponential"){
    stopifnot(all(nat_mat < 0))
    -1/nat_mat
  } else if(family == "curved_gaussian"){
    stopifnot(all(nat_mat > 0))
    1/nat_mat
  } else {
    stop("family not found")
  }
}

######

.compute_mean_neg_binom <- function(nat_mat, scalar, ...){
  if(is.na(scalar)) stop("No argument scalar provided for negative binomial")
  stopifnot(length(scalar) == 1)

  scalar*exp(nat_mat)/(1-exp(nat_mat))
}

.dictate_direction <- function(family){
  direction <- NA

  if(family %in% c("exponential", "neg_binom")) {
    direction <- "<="
  } else if(family %in% c("curved_gaussian")) {
    direction <- ">="
  } else if(!family %in% c("gaussian", "poisson")) {
    stop("family not found")
  }

  direction
}

.mean_transformation <- function(dat, family, tol = 1e-3, ...){
  if(family == "exponential"){
    dat <- -1/(dat+1)
  } else if(family == "curved_gaussian"){
    dat <- 1/(dat+1)
  } else if(family == "poisson"){
    dat <- log(dat + tol)
  } else if(family == "neg_binom"){
    dat <- .mean_transformation_neg_binom(dat, tol, ...)
  } else if(family != "gaussian") {
    stop("family not found")
  }

  dat
}

.mean_transformation_neg_binom <- function(dat, tol, scalar = NA, ...){
  if(is.na(scalar)) stop("No argument scalar provided for negative binomial")
  stopifnot(length(scalar) == 1)

  dat_new <- (dat + tol)/scalar
  log(dat_new / (1+dat_new))
}

.compute_prediction_interval_from_mean <- function(val, family, width, ...){
  stopifnot(width >= 0, width <= 1)

  lower_val <- 0.5-width/2; upper_val <- 0.5+width/2

  if(family == "exponential"){
    c(stats::qexp(lower_val, rate = 1/val), stats::qexp(upper_val, rate = 1/val))
  } else if(family == "poisson"){
    c(stats::qpois(lower_val, lambda = val), stats::qpois(upper_val, lambda = val))
  } else if(family == "gaussian"){
    .compute_prediction_interval_gaussian(val, width, lower_val, upper_val, ...)
  } else if(family == "neg_binom"){
    .compute_prediction_interval_neg_binom(val, width, lower_val, upper_val,  ...)
  } else {
    .compute_prediction_interval_curved_gaussian(val, width, lower_val, upper_val, ...)
  }
}

.compute_prediction_interval_gaussian <- function(val, width, lower_val, upper_val, scalar, ...){
  if(is.na(scalar)) stop("No argument scalar provided for gaussian")

  c(stats::qnorm(lower_val, mean = val, sd = scalar), stats::qnorm(upper_val, mean = val, sd = scalar))
}

.compute_prediction_interval_neg_binom <- function(val, width, lower_val, upper_val, scalar, ...){
  if(is.na(scalar)) stop("No argument scalar provided for negative binomial")

  c(stats::qnbinom(lower_val, size = scalar, mu = val), stats::qnbinom(upper_val, size = scalar, mu = val))
}

.compute_prediction_interval_curved_gaussian <- function(val, width, lower_val, upper_val, scalar, ...){
  if(is.na(scalar)) stop("No argument scalar provided for curved gaussian")

  c(stats::qnorm(lower_val, mean = val, sd = val/scalar), stats::qnorm(upper_val, mean = val, sd = val/scalar))
}
