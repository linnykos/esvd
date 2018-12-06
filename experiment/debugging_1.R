rm(list=ls())
load("../results/tmp.RData")
i=182
.em_mixture(dat[,i], mixture = "gamma.tgaussian", min_val = 0)

############

x = dat[,i]
mixture = "gamma.tgaussian"
min_val = 0
max_iter = 100

param <- .initialize_mixture(x, mixture, min_val)
max_prop <- param[["proportion"]]

x2 <- .jitter_zeros(x)
eps <- 10
iter <- 0
loglik_old <- 0
#
# while(eps > 0.5) {
#   print(iter)
#   #E-step
#   wt <- .calculate_weight(x2, param)
#
#   #M-step
#   param["proportion"] <- min(sum(wt[, 1])/nrow(wt), max_prop)
#   param[["class1"]] <- estimate_parameter(param[["class1"]], x2, wt[,1])
#   param[["class2"]] <- estimate_parameter(param[["class2"]], x, wt[,2],
#                                           min_mean = max(0, compute_mean(param[["class1"]])))
#
#   #see if converged
#   loglik <- .log_likelihood(x2, param)
#   eps <- (loglik - loglik_old)^2
#   loglik_old <- loglik
#   iter <- iter + 1
#   if (iter > max_iter) break()
# }

############

print(iter)
#E-step
wt <- .calculate_weight(x2, param)

#M-step
param["proportion"] <- min(sum(wt[, 1])/nrow(wt), max_prop)
param[["class1"]] <- estimate_parameter(param[["class1"]], x2, wt[,1])

