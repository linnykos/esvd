likelihood <- function(vec, x){
  UseMethod("likelihood")
}

likelihood.exponential <- function(vec, x, min_val = 0){
  stats::dexp(x - min_val, rate = vec["rate"])
}

likelihood.gamma <- function(vec, x){
  stats::dgamma(x, shape = vec["shape"], rate = vec["rate"])
}

likelihood.gaussian <- function(vec, x){
  stats::dnorm(x, mean = vec["mean"], sd = vec["sd"])
}

likelihood.tgaussian  <- function(vec, x){
  likelihood.gaussian(vec, x)/(1 - stats::pnorm(0, mean = vec["mean"], sd = vec["sd"]))
}

##########

compute_mean <- function(obj){
  UseMethod("compute_mean")
}

compute_mean.exponential <- function(obj){
  1/obj["rate"]
}

compute_mean.gamma <- function(obj){
  obj["shape"]/obj["rate"]
}

compute_mean.gaussian <- function(obj){
  obj["mean"]
}

compute_mean.tgaussian <- function(obj){
  obj["mean"]
}

#########

initialize.exponential <- function(x){
  structure(c(rate = 200), class = "exponential")
}

initialize.gamma <- function(x){
  structure(c(shape = 1, rate = 200), class = "gamma")
}

initialize.gaussian <- function(x, min_val = 0){
  structure(c(mean = mean(x[x > min_val]), sd = stats::sd(x[x > min_val])), class = "gaussian")
}

initialize.tgaussian <- function(x, min_val = 0){
  structure(c(mean = mean(x[x > min_val]), sd = stats::sd(x[x > min_val])), class = "tgaussian")
}

#########

estimate_parameter <- function(obj, x, weight = rep(1, length(x)), ...){
  UseMethod("estimate_parameter")
}

estimate_parameter.exponential <- function(obj, x, weight = rep(1, length(x)),
                                           min_val = 0, ...){
  mean_val <- sum(weight * (x - min_val))/sum(weight)
  structure(c(rate = 1/mean_val), class = "exponential")
}

estimate_parameter.gamma <- function(obj, x, weight = rep(1, length(x)), ...) {
  fn <- function(alpha, target){
    log(alpha) - digamma(alpha) - target
  }

  tp_s <- sum(weight)
  tp_t <- sum(weight * x)
  tp_u <- sum(weight * log(x))
  tp_v <- -tp_u / tp_s - log(tp_s / tp_t)
  if (tp_v <= 0){
    alpha <- 20
  }else{
    alpha0 <- (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
    if (alpha0 >= 20){alpha = 20
    }else{
      alpha <- stats::uniroot(fn, c(0.9, 1.1) * alpha0, target = tp_v,
                              extendInt = "yes")$root
    }
  }
  beta <- tp_s / tp_t * alpha
  structure(c(shape = alpha, rate = beta), class = "gamma")
}

estimate_parameter.gaussian <- function(obj, x, weight = rep(1, length(x)), ...){
 mean_val <- sum(weight * x)/sum(weight)
 sd_val <- sqrt(sum(weight * (x - mean_val)^2)/sum(weight))

 structure(c(mean = mean_val, sd = sd_val), class = "gaussian")
}


estimate_parameter.tgaussian <- function(obj, x, weight = rep(1, length(x)),
                                         min_val = 0,
                                         min_mean = 0){
  func <- function(vec){
    -(sum(weight * stats::dnorm(x, vec[1], vec[2], log = TRUE)) -
        length(x) * log(1-stats::pnorm(0, mean = vec[1], sd = vec[2])))
  }

  if(all(x == min_val)) return(c(mean = NA, sd = NA))

  min_sd <- 2*stats::sd(x[x > min_val])

  res <- stats::optim(par = c(mean(x[x > min_val]), 2*stats::sd(x[x > min_val])),
               fn = func, method = "L-BFGS-B", lower = c(min_mean, min_sd))$par

  structure(c(mean = res[1], sd = res[2]), class = "tgaussian")
}
