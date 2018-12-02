rm(list=ls())
load("../results/wasserstein_simulation.RData")

#forbenius_loss <- sapply(1:length(res), function(i){
loss_mat <- sapply(1:length(res), function(i){
  print(i)
  sapply(1:trials, function(y){
    if(length(res[[i]][[y]]) == 1) return(rep(NA,2))
    cat('*')
    set.seed(y)
    vec <- paramMat[i,]
    obj <- .data_generator(cell_pop, gene_pop, n_each = vec["n"], d_each = vec["d"])
    n <- nrow(obj$cell_mat)
    forbenius <- .l2norm(obj$cell_mat - res[[i]][[y]]$u_mat)^2/n
    wasserstein <- transport::wasserstein(transport::pp(obj$cell_mat),
                                          transport::pp(res[[i]][[y]]$u_mat), p = 1)

    c(forbenius, wasserstein)
  })
})

save.image("../results/wasserstein_simulation_values.RData")

forbenius_vec <- sapply(loss_mat, function(x){
  median(x[1,], na.rm = T)
})
bound_vec <- apply(forbenius_loss, 2, stats::quantile, na.rm = T,
                   probs = c(0.25,0.75))
plot(paramMat[,"n"], forbenius_vec, ylim = range(bound_vec))
points(paramMat[,"n"], bound_vec[1,], col = "red", pch = 16)
points(paramMat[,"n"], bound_vec[2,], col = "blue", pch = 16)
