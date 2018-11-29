rm(list=ls())
load("../results/wasserstein_simulation.RData")

forbenius_loss <- sapply(1:length(res), function(i){
  print(i)
  sapply(1:trials, function(y){
    cat('*')
    set.seed(y)
    vec <- paramMat[i,]
    obj <- .data_generator(cell_pop, gene_pop, n_each = vec["n"], d_each = vec["d"])
    n <- nrow(obj$cell_mat)
    .l2norm(obj$cell_mat - res[[i]][[y]]$u_mat)/n
  })
})

forbenius_vec <- apply(forbenius_loss, 2, mean)
