rm(list=ls())
load("../results/factorization_results.RData")

ncores <- 15
doMC::registerDoMC(cores = ncores)

func <- function(i){
  print(i)

  sapply(1:6, function(k){
    transport::wasserstein(transport::pp(res[[1]][[i]]$dat$truth),
                           transport::pp(res[[1]][[i]][[(k-1)*2+1]]), p = 1)
  })
}

res_wasserstein <- foreach::"%dopar%"(foreach::foreach(i = 1:trials), func(i))

save.image("../results/factorization_results_wasserstein.RData")
