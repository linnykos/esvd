source("../experiment/em_gamma_normal.R")
source("../experiment/em_gamma_truncatednormal.R")

### estimate parameters in the mixture distribution
get_mix_switch = function(xdata, point = log10(1.01), prop_init = NA,
                          truncated = F){
  if(!truncated){
    paramt <- get_mix(xdata, point = point, prop_init = prop_init)
  } else {
    paramt <- .get_mix(xdata, point = point, prop_init = prop_init)
  }

  paramt <- structure(c("proportion" = paramt[1],
                        "shape" = paramt[2],
                        "rate" = paramt[3],
                        "mean" = paramt[4],
                        "sd" = paramt[5]))
  if(!truncated) class(paramt) <- "Gamma-Normal" else class(paramt) <- "Gamma-TNormal"


  xdata_rm <- xdata[xdata > point]

  #evaluate likelihood of only normal component
  lik <- sum(sapply(xdata_rm, function(x){
    val1 <- (paramt[1])*dgamma(x, shape = paramt[2], rate = paramt[3])
    val2 <- (1-paramt[1])*dnorm(x, mean = paramt[4], sd = paramt[5])
    if(truncated) val2 <- val2 / (1 - stats::pnorm(0, mean = paramt[4], sd = paramt[5]))
    val1 + val2
  }))

  #compare
  xdata_rm <- xdata[xdata > point]
  mu <- mean(xdata_rm); sig <- sd(xdata_rm)
  prop <- (length(xdata) - length(xdata_rm))/length(xdata)
  lik2 <- sum(sapply(xdata, function(x){
    (1-prop)*dnorm(x, mean = mu, sd = sig)
  }))

  if(lik2 > lik){
    paramt <- structure(c("proportion" = prop,
                         "delta" = point,
                         "mean" = mu,
                         "sd" = sig), class = "Delta-Normal")
  }

  return(paramt)
}
