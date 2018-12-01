rm(list=ls())
load("../results/wasserstein_simulation.RData")

#forbenius_loss <- sapply(1:length(res), function(i){
forbenius_loss <- sapply(1:length(res), function(i){
  print(i)
  sapply(1:trials, function(y){
    if(length(res[[i]][[y]]) == 1) return(NA)
    cat('*')
    set.seed(y)
    vec <- paramMat[i,]
    obj <- .data_generator(cell_pop, gene_pop, n_each = vec["n"], d_each = vec["d"])
    n <- nrow(obj$cell_mat)
    .l2norm(obj$cell_mat - res[[i]][[y]]$u_mat)^2/n
  })
})

forbenius_vec <- apply(forbenius_loss, 2, median, na.rm = T)
bound_vec <- apply(forbenius_loss, 2, stats::quantile, na.rm = T,
                   probs = c(0.25,0.75))
plot(paramMat[,"n"], forbenius_vec, ylim = range(bound_vec))
points(paramMat[,"n"], bound_vec[1,], col = "red", pch = 16)
points(paramMat[,"n"], bound_vec[2,], col = "blue", pch = 16)
