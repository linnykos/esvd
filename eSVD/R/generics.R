.evaluate_objective <- function (dat, u_mat, v_mat, ...) {
  stopifnot(is.matrix(dat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat),
            ncol(u_mat) == ncol(v_mat))

  if(attr(dat, "family") == "gaussian"){
    .evaluate_objective.gaussian(dat, u_mat, v_mat, ...)
  } else if(attr(dat, "family") == "poisson"){
    .evaluate_objective.poisson(dat, u_mat, v_mat, ...)
  } else if(attr(dat, "family") == "neg_binom"){
    .evaluate_objective.neg_binom(dat, u_mat, v_mat, ...)
  } else if(attr(dat, "family") == "exponential"){
    .evaluate_objective.exponential(dat, u_mat, v_mat, ...)
  } else if(attr(dat, "family") == "curved_gaussian"){
    .evaluate_objective.curved_gaussian(dat, u_mat, v_mat, ...)
  } else {
    stop()
  }
}

##########

.evaluate_objective_single <- function (dat_vec, current_vec, other_mat, n, p, ...) {
  stopifnot(!is.matrix(dat_vec))
  stopifnot(length(current_vec) == ncol(other_mat))
  stopifnot(length(dat_vec) == nrow(other_mat))

  if(attr(dat_vec, "family") == "gaussian"){
    .evaluate_objective_single.gaussian(dat_vec, current_vec, other_mat, n, p, ...)
  } else if(attr(dat_vec, "family") == "poisson"){
    .evaluate_objective_single.poisson(dat_vec, current_vec, other_mat, n, p, ...)
  } else if(attr(dat_vec, "family") == "neg_binom"){
    .evaluate_objective_single.neg_binom(dat_vec, current_vec, other_mat, n, p, ...)
  } else if(attr(dat_vec, "family") == "exponential"){
    .evaluate_objective_single.exponential(dat_vec, current_vec, other_mat, n, p, ...)
  } else if(attr(dat_vec, "family") == "curved_gaussian"){
    .evaluate_objective_single.curved_gaussian(dat_vec, current_vec, other_mat, n, p, ...)
  } else {
    stop()
  }
}

###########

.gradient_vec <- function (dat_vec, current_vec, other_mat, ...) {
  stopifnot(!is.matrix(dat_vec))
  stopifnot(length(current_vec) == ncol(other_mat))
  stopifnot(length(dat_vec) == nrow(other_mat))

  if(attr(dat_vec, "family") == "gaussian"){
    .gradient_vec.gaussian(dat_vec, current_vec, other_mat, ...)
  } else if(attr(dat_vec, "family") == "poisson"){
    .gradient_vec.poisson(dat_vec, current_vec, other_mat, ...)
  } else if(attr(dat_vec, "family") == "neg_binom"){
    .gradient_vec.neg_binom(dat_vec, current_vec, other_mat, ...)
  } else if(attr(dat_vec, "family") == "exponential"){
    .gradient_vec.exponential(dat_vec, current_vec, other_mat, ...)
  } else if(attr(dat_vec, "family") == "curved_gaussian"){
    .gradient_vec.curved_gaussian(dat_vec, current_vec, other_mat, ...)
  } else {
    stop()
  }
}

###########

.evaluate_objective_mat <- function(dat, nat_mat, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)))

  if(attr(dat, "family") == "gaussian"){
    .evaluate_objective_mat.gaussian(dat, nat_mat, ...)
  } else if(attr(dat, "family") == "poisson"){
    .evaluate_objective_mat.poisson(dat, nat_mat, ...)
  } else if(attr(dat, "family") == "neg_binom"){
    .evaluate_objective_mat.neg_binom(dat, nat_mat, ...)
  } else if(attr(dat, "family") == "exponential"){
    .evaluate_objective_mat.exponential(dat, nat_mat, ...)
  } else if(attr(dat, "family") == "curved_gaussian"){
    .evaluate_objective_mat.curved_gaussian(dat, nat_mat, ...)
  } else {
    stop()
  }
}

###########

#' Gradient of the objective function
#'
#' Computes the gradient for a particular model (based on
#' \code{attr(dat, "family")}) of the \code{.evaluate_objective_mat} function.
#'
#' Note, \code{dat} is NOT allowed to have any \code{NA} values for this
#' function.
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param nat_mat \code{n} by \code{d} matrix where each entry represents the natural parameter of the corresponding entry in \code{dat}
#' @param ... other parameters
#'
#' @return \code{n} by \code{d} matrix
.gradient_mat <- function(dat, nat_mat, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)))

  if(attr(dat, "family") == "gaussian"){
    .gradient_mat.gaussian(dat, nat_mat, ...)
  } else if(attr(dat, "family") == "poisson"){
    .gradient_mat.poisson(dat, nat_mat, ...)
  } else if(attr(dat, "family") == "neg_binom"){
    .gradient_mat.neg_binom(dat, nat_mat, ...)
  } else if(attr(dat, "family") == "exponential"){
    .gradient_mat.exponential(dat, nat_mat, ...)
  } else if(attr(dat, "family") == "curved_gaussian"){
    .gradient_mat.curved_gaussian(dat, nat_mat, ...)
  } else {
    stop()
  }
}
