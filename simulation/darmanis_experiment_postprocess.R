rm(list=ls())
load("../experiment/darmanis_experiment.RData")

len <- nrow(paramMat)
for(i in 1:len){
  plot_mat <- sapply(1:trials, function(j){
    cbind(res[[i]][[j]]$missing_val, res[[i]][[j]]$pred_val)
  })

  plot_mat <- do.call(rbind, plot_mat)

  plot(plot_mat[,1], plot_mat[,2], asp = T, main = i, pch = 16)
  lines(c(-1e3,1e3), c(-1e3, 1e3), col = "red", lty = 2, lwd = 2)
}
