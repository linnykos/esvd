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
    dat <- -1/dat
  } else if(family == "curved_gaussian"){
    dat <- 1/dat
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
  if(is.na(size)) stop("No argument r provided for negative binomial")

  dat_new <- (dat + tol)/size
  log(dat_new / (1+dat_new))
}
