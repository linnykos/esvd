rm(list=ls())

load("../results/factorization_esvd.RData")

# plot the embeddings
k <- 1
par(mfrow = c(3,3))
for(i in 1:8){
  plot(res[[(k-1)*8+i]][[1]]$fit[,1], res[[(k-1)*8+i]][[1]]$fit[,2], asp = T, col = rep(1:4, each = paramMat[1,"n_each"]),
       pch = 16, main = i)
}
