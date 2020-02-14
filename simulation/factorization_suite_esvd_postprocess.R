rm(list=ls())

load("../results/factorization_esvd.RData")

# plot the embeddings
k <- 1
par(mfrow = c(3,3))
for(i in 1:8){
  plot(res[[(k-1)*8+i]][[1]]$fit[,1], res[[(k-1)*8+i]][[1]]$fit[,2], asp = T, col = rep(1:4, each = paramMat[1,"n_each"]),
       pch = 16, main = i)
}

# plot the missing values
k <- 1
par(mfrow = c(3,3))
for(i in 1:8){
  res[[(k-1)*8+i]] <- res[[(k-1)*8+i]][which(sapply(res[[(k-1)*8+i]], function(x){!all(is.na(x))}))]

  plot_mat <- lapply(1:length(res[[(k-1)*8+i]]), function(j){
    cbind(res[[(k-1)*8+i]][[j]]$missing_val, res[[(k-1)*8+i]][[j]]$pred_val)
  })

  if(length(plot_mat) > 1) plot_mat <- do.call(rbind, plot_mat) else plot_mat <- plot_mat[[1]]

  plot(plot_mat[,1], plot_mat[,2], asp = T, main = i, pch = 16)
  lines(c(-1e5,1e5), c(-1e5, 1e5), col = "red", lty = 2, lwd = 2)
}
