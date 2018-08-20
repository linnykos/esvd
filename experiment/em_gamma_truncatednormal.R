### root-finding equation
.fn = function(alpha, target){
  log(alpha) - digamma(alpha) - target
}

### update parameters in gamma distribution
.update_gmm_pars = function(x, wt) {
  tp_s = sum(wt)
  tp_t = sum(wt * x)
  tp_u = sum(wt * log(x))
  tp_v = -tp_u / tp_s - log(tp_s / tp_t)
  if (tp_v <= 0){
    alpha = 20
  }else{
    alpha0 = (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v)) / 12 / tp_v
    if (alpha0 >= 20){alpha = 20
    }else{
      alpha = stats::uniroot(.fn, c(0.9, 1.1) * alpha0, target = tp_v,
                             extendInt = "yes")$root
    }
  }
  ## need to solve log(x) - digamma(x) = tp_v
  ## We use this approximation to compute the initial value
  beta = tp_s / tp_t * alpha
  return(c(alpha, beta))
}

### estimate parameters in the mixture distribution
.get_mix = function(xdata, point = log10(1.01)){
  #initialize the five parameters
  inits = rep(0, 5)
  inits[1] = sum(xdata == point)/length(xdata)
  if (inits[1] == 0) {inits[1] = 0.01}
  inits[2:3] = c(0.5, 1)
  xdata_rm = xdata[xdata > point]
  inits[4:5] = c(mean(xdata_rm), sd(xdata_rm))
  if (is.na(inits[5])) {inits[5] = 0}
  paramt = inits
  eps = 10
  iter = 0
  loglik_old = 0

  while(eps > 0.5) {
    #E-step
    wt = .calculate_weight(xdata, paramt)

    #M-step
    paramt[1] = sum(wt[, 1])/nrow(wt)
    res <- estimate_truncated_normal(xdata, wt[,2])
    paramt[4] = res$mean
    paramt[5] = sqrt(res$var)
    paramt[2:3] = .update_gmm_pars(x=xdata, wt=wt[,1])

    #see if converged
    loglik = sum(log10(.dmix(xdata, paramt)))
    eps = (loglik - loglik_old)^2
    loglik_old = loglik
    iter = iter + 1
    if (iter > 100)
      break
  }
  return(paramt)
}

.dmix <- function (x, pars) {
  pars[1] * stats::dgamma(x, shape = pars[2], rate = pars[3]) +
    (1 - pars[1]) * stats::dnorm(x, mean = pars[4], sd = pars[5])/(1-stats::pnorm(0))
}

.calculate_weight <- function (x, paramt) {
  pz1 = paramt[1] * stats::dgamma(x, shape = paramt[2], rate = paramt[3])
  pz2 = (1 - paramt[1]) * stats::dnorm(x, mean = paramt[4], sd = paramt[5])/(1-stats::pnorm(0))
  pz = pz1/(pz1 + pz2)
  pz[pz1 == 0] = 0
  return(cbind(pz, 1 - pz))
}
