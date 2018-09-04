.calculate_weight <- function (x, param_list) {
  pz1 <- param_list[["proportion"]] * likelihood(param_list[["class1"]], x)
  pz2 <- (1 - param_list[["proportion"]]) * likelihood(param_list[["class2"]], x)
  pz <- pz1/(pz1 + pz2)
  pz[pz1 == 0] <- 0

  cbind(pz, 1 - pz)
}

.log_likelihood <- function(x, param_list){
  sum(log10(param_list[["proportion"]] * likelihood(param_list[["class1"]], x) +
              (1-param_list[["proportion"]]) * likelihood(param_list[["class2"]], x)))
}

.em_mixture = function(x, mixture = "gamma.tgaussian",
                       min_val = log10(1.01), prop_init = NA, max_iter = 100){
  param_list <- .initialize_mixture(x, mixture, min_val, prop_init)

  x2 <- .jitter_zeros(x)
  eps <- 10
  iter <- 0
  loglik_old <- 0

  while(eps > 0.5) {
    #E-step
    wt <- .calculate_weight(x, param_list)

    #M-step
    param_list["proportion"] <- max(sum(wt[, 1])/nrow(wt), 0.1)
    param_list[["class1"]] <- estimate_parameter(param_list[["class1"]], x2, wt[,1])
    param_list[["class2"]] <- estimate_parameter(param_list[["class2"]], x, wt[,2],
                                                 min_mean = compute_mean(param_list[["class1"]]))

    #see if converged
    loglik <- .log_likelihood(x, param_list)
    eps <- (loglik - loglik_old)^2
    loglik_old <- loglik
    iter <- iter + 1
    if (iter > max_iter) break()
  }

  param_list
}

.initialize_mixture <- function(x, mixture, min_val, prop_init){
  proportion <- length(which(x == min_val))/length(x)
  class_vec <- strsplit(mixture, split = "\\.")[[1]]
  lis <- lapply(class_vec, function(y){
    eval(parse(text = paste0("initialize.", y, "(x)")))
  })

  structure(list(class1 = lis[[1]], class2 = lis[[2]], proportion = proportion))
}

.jitter_zeros <- function(x, min_val = log10(1.01)){
  idx <- which(x == min_val)
  if(length(idx) == 0) return(x)

  dif <- min(x[x > min_val]) - min_val
  x[idx] <- x[idx] + stats::rexp(length(idx), rate = 2/dif)

  x
}


