.dictate_direction <- function(family = "exponential"){
  if(family %in% c("exponential")){
    direction <- "<="
  } else if(family %in% c("gaussian", "poisson")){
    direction <- ">="
  } else {
    stop("family not found")
  }

  direction
}

.mean_transformation <- function(dat, family = "exponential", tol = 1e-3,
                                 scalar = 2){
  if(family == "exponential"){
    dat <- -1/dat
  } else if(family == "gaussian"){
    dat <- scalar^2/dat
  } else if(family == "poisson"){
    dat <- log(dat + tol)
  } else {
    stop("family not found")
  }

  dat
}
