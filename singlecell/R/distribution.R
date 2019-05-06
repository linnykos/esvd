.dictate_direction <- function(family = "exponential"){
  if(family %in% c("exponential", "gaussian")) {
    direction <- "<="
  } else if(family %in% c("poisson")) {
    direction <- ">="
  } else {
    stop("family not found")
  }

  direction
}

.mean_transformation <- function(dat, family = "exponential", tol = 1e-3){
  if(family == "exponential"){
    dat <- -1/dat
  } else if(family == "gaussian"){
    dat <- -1/dat
  } else if(family == "poisson"){
    dat <- log(dat + tol)
  } else {
    stop("family not found")
  }

  dat
}
