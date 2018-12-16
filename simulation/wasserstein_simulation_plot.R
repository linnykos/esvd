rm(list=ls())
load("../results/wasserstein_simulation.RData")

#####################

# zz <- res[[10]][[1]]$u_mat %*% t(res[[10]][[1]]$v_mat)
# tmpu <- svd(zz); tmpu <- (nrow(res[[10]][[1]]$u_mat)/nrow(res[[10]][[1]]$v_mat))^(1/4)*tmpu$u %*% diag(sqrt(tmpu$d))
# eigen_val <- eigen(stats::cov(tmpu))$values
# #only 2 eigenvalues
# eigen_val <- eigen_val[1:2]
# range(log(eigen_val)/(-log(c(1.01,2))))

#forbenius_loss <- sapply(1:length(res), function(i){
loss_mat <- lapply(1:length(res), function(i){
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

forbenius_bound <- sapply(loss_mat, function(x){
  stats::quantile(x[1,], na.rm = T, probs = c(0.25, 0.5, 0.75))
})

wasserstein_bound <- sapply(loss_mat, function(x){
  stats::quantile(x[2,], na.rm = T,probs = c(0.25, 0.5, 0.75))
})

#########################

# figure out hypothetical bounds
vec1 <- forbenius_bound[2,]
vec2 <- 4*paramMat[,"n"]
range(-log(vec1)/log(vec2))

png("../figure/simulation/wasserstein.png", height = 1100, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,4,4,0.5))
plot(paramMat[,"n"]*4, forbenius_bound[2,], ylim = range(forbenius_bound),
     pch = 16, xlab = "n", ylab = "Forbenius loss",
     main = "Forbenius loss for cells")
lines(paramMat[,"n"]*4, forbenius_bound[2,], lwd = 2)
Hmisc::errbar(paramMat[,"n"]*4, forbenius_bound[2,], yplus = forbenius_bound[3,],
              yminus = forbenius_bound[1,], add = T)

plot(paramMat[,"n"]*4, wasserstein_bound[2,], ylim = range(wasserstein_bound),
     pch = 16, xlab = "n", ylab = "Wasserstein loss",
     main = "Wasserstein loss for cells")
lines(paramMat[,"n"]*4, wasserstein_bound[2,], lwd = 2)
Hmisc::errbar(paramMat[,"n"]*4, wasserstein_bound[2,], yplus = wasserstein_bound[3,],
              yminus = wasserstein_bound[1,], add = T)
graphics.off()
