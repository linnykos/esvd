rm(list=ls())
load("../results/factorization_results.RData")

ncores <- 15
doMC::registerDoMC(cores = ncores)

func <- function(i){
  print(i)

  sapply(1:6, function(k){
    # try all orientations of the data
    dat_truth <- res[[1]][[i]]$dat$truth[,c(2,1)]
    dat_truth <- apply(dat_truth, 2, function(x){(x-mean(x))/sd(x)})
    dat_est <- res[[1]][[i]][[(k-1)*2+1]]
    dat_est <- apply(dat_est, 2, function(x){(x-mean(x))/sd(x)})
    est_1 <- dat_est[,1]
    est_2 <- dat_est[,2]

    val1 <- transport::wasserstein(transport::pp(dat_truth),
                                   transport::pp(cbind(est_1, est_2)), p = 1)
    val2 <- transport::wasserstein(transport::pp(dat_truth),
                                   transport::pp(cbind(-est_1, est_2)), p = 1)
    val3 <- transport::wasserstein(transport::pp(dat_truth),
                                   transport::pp(cbind(est_1, -est_2)), p = 1)
    val4 <- transport::wasserstein(transport::pp(dat_truth),
                                   transport::pp(cbind(-est_1, -est_2)), p = 1)

    min(c(val1, val2, val3, val4))
  })
}

res_wasserstein <- foreach::"%dopar%"(foreach::foreach(i = 1:trials), func(i))

save.image("../results/factorization_results_wasserstein.RData")
