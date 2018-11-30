rm(list=ls())
load("../results/wasserstein_simulation.RData")
load("../results/wasserstein_tmp.RData")

#forbenius_loss <- sapply(1:length(res), function(i){
forbenius_loss <- sapply(1:8, function(i){
  print(i)
#  sapply(1:trials, function(y){
  sapply(1:50, function(y){
    if(length(res[[i]][[y]]) == 1) return(NA)
    cat('*')
    set.seed(y)
    vec <- paramMat[i,]
    obj <- .data_generator(cell_pop, gene_pop, n_each = vec["n"], d_each = vec["d"])
    n <- nrow(obj$cell_mat)
    .l2norm(obj$cell_mat - res[[i]][[y]]$u_mat)^2/n
  })
})

forbenius_vec <- apply(forbenius_loss, 2, mean, na.rm = T)
