rm(list=ls())
load("../results/wasserstein_simulation.RData")

######### # test
tmp <- res[[7]][[1]]$res_our
plot(tmp[,1], tmp[,2], asp = T, col = rep(1:4, each = nrow(tmp)/4), pch = 16)

#####################

# collect the values
loss_mat <- matrix(0, nrow = 2*nrow(paramMat), ncol = trials)
for(i in 1:nrow(paramMat)){
  for(j in 1:trials){
    loss_mat[(i-1)+1,j] <- res[[i]][[j]]$l2_loss
    loss_mat[2*i,j] <- res[[i]][[j]]$wasserstein_loss
  }
}

l2_val <- matrix(0, nrow = nrow(paramMat), ncol = 3)
for(i in 1:nrow(paramMat)){
  #mean_val <- median(loss_mat[(i-1)+1,])
  #sd_val <- sd(loss_mat[(i-1)+1,])
  #l2_val[i,] <- c(mean_val-sd_val, mean_val, mean_val+sd_val)
  l2_val[i,] <- c(quantile(loss_mat[(i-1)+1,], probs = 0.25),
                  quantile(loss_mat[(i-1)+1,], probs = 0.5),
                  quantile(loss_mat[(i-1)+1,], probs = 0.75))
}

wass_val <- matrix(0, nrow = nrow(paramMat), ncol = 3)
for(i in 1:nrow(paramMat)){
  mean_val <- median(loss_mat[2*i,])
  sd_val <- sd(loss_mat[2*i,])
  wass_val[i,] <- c(mean_val-sd_val, mean_val, mean_val+sd_val)
}

#########################

keep_idx <- 1:9
png("../figure/simulation/wasserstein.png", height = 1200, width = 1500, res = 300, units = "px")
plot(paramMat[keep_idx,"n_each"]*4, l2_val[keep_idx,2], ylim = range(l2_val[keep_idx,]),
     pch = 16, xlab = "n", ylab = "Forbenius loss",
     main = "Forbenius loss (squared) for cells")
lines(paramMat[keep_idx,"n_each"]*4, l2_val[keep_idx,2], lwd = 2)
Hmisc::errbar(paramMat[keep_idx,"n_each"]*4, l2_val[keep_idx,2], yplus = l2_val[keep_idx,3],
              yminus = l2_val[keep_idx,1], add = T)
graphics.off()
