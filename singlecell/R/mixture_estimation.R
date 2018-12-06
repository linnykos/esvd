#' EM mixture for modeling dropout
#'
#' @param x vector
#' @param mixture family, represented by a string
#' @param min_val minimum value used for initialization
#' @param max_iter maximum iteration for EM algorithm
#'
#' @return parameter object
#' @export
em_mixture <- function(x, mixture = "gamma.tgaussian",
                        min_val = 0, max_iter = 100){
  param <- .initialize_mixture(x, mixture, min_val)
  max_prop <- param[["proportion"]]

  x2 <- .jitter_zeros(x)
  eps <- 10
  iter <- 0
  loglik_old <- 0

  while(eps > 0.5) {
    #E-step
    wt <- .calculate_weight(x2, param)

    #M-step
    param["proportion"] <- min(sum(wt[, 1])/nrow(wt), max_prop)
    param[["class1"]] <- estimate_parameter(param[["class1"]], x2, wt[,1])
    param[["class2"]] <- estimate_parameter(param[["class2"]], x, wt[,2],
                                            min_mean = max(0, compute_mean(param[["class1"]])))

    #see if converged
    loglik <- .log_likelihood(x2, param)
    eps <- (loglik - loglik_old)^2
    loglik_old <- loglik
    iter <- iter + 1
    if (iter > max_iter) break()
  }

  param
}

.calculate_weight <- function (x, param) {
  pz1 <- param[["proportion"]] * likelihood(param[["class1"]], x)
  pz2 <- (1 - param[["proportion"]]) * likelihood(param[["class2"]], x)
  pz <- pz1/(pz1 + pz2)
  pz[pz1 == 0] <- 0

  cbind(pz, 1 - pz)
}

.log_likelihood <- function(x, param){
  sum(log10(param[["proportion"]] * likelihood(param[["class1"]], x) +
              (1-param[["proportion"]]) * likelihood(param[["class2"]], x)))
}

.initialize_mixture <- function(x, mixture, min_val){
  proportion <- max(length(which(x == min_val))/length(x), 0.05)
  class_vec <- strsplit(mixture, split = "\\.")[[1]]
  lis <- lapply(class_vec, function(y){
    eval(parse(text = paste0("initialize.", y, "(x)")))
  })

  structure(list(class1 = lis[[1]], class2 = lis[[2]], proportion = proportion))
}

.jitter_zeros <- function(x, min_val = 0){
  idx <- which(x <= min_val)
  if(length(idx) == 0) return(x)
  if(length(idx) == length(x)) return(x)

  dif <- min(x[-idx]) - min_val
  x[idx] <- x[idx] + stats::rexp(length(idx), rate = 2/dif)

  stopifnot(length(which(x == min_val)) == 0)

  x
}


