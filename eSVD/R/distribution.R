#' Compute mean from natural parameter matrix
#'
#' @param nat_mat matrix
#' @param family character
#' @param ... additional parameters for distribution
#'
#' @return matrix
#' @export
compute_mean <- function(nat_mat, family, ...){
  if(family == "gaussian") {
    return(nat_mat)
  } else if(family == "poisson"){
    exp(nat_mat)
  } else if(family == "neg_binom"){
    .compute_mean_neg_binom(nat_mat, ...)
  } else if(family == "exponential"){
    -1/nat_mat
  } else if(family == "curved_gaussian"){
    1/nat_mat
  } else {
    stop("family not found")
  }
}

######

.compute_mean_neg_binom <- function(nat_mat, size){
  if(is.na(size)) stop("No argument size provided for negative binomial")

  size*exp(nat_mat)/(1-exp(nat_mat))
}

.dictate_direction <- function(family){
  direction <- NA

  if(family %in% c("exponential", "neg_binom")) {
    direction <- "<="
  } else if(family %in% c("poisson", "curved_gaussian")) {
    direction <- ">="
  } else if(family != "gaussian") {
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

.mean_transformation_neg_binom <- function(dat, tol, size = NA){
  if(is.na(size)) stop("No argument size provided for negative binomial")

  dat_new <- (dat + tol)/size
  log(dat_new / (1+dat_new))
}


