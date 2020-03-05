rm(list=ls())

load("../results/factorization_esvd.RData")

k = 4
i = 9
plot_mat <- lapply(1:length(res[[(k-1)*8+i]]), function(j){
  cbind(res[[(k-1)*8+i]][[j]]$missing_val, res[[(k-1)*8+i]][[j]]$pred_val)
})

if(length(plot_mat) > 1) plot_mat <- do.call(rbind, plot_mat) else plot_mat <- plot_mat[[1]]

(nrow(plot_mat)/(sum((plot_mat[,1]-plot_mat[,2])^2/(plot_mat[,2])^2)))^(1/2)

##################

k = 3
i = 2
plot_mat <- lapply(1:length(res[[(k-1)*8+i]]), function(j){
  cbind(res[[(k-1)*8+i]][[j]]$missing_val, res[[(k-1)*8+i]][[j]]$pred_val)
})

if(length(plot_mat) > 1) plot_mat <- do.call(rbind, plot_mat) else plot_mat <- plot_mat[[1]]

r_seq <- 1:200
vec <- sapply(r_seq, function(x){
  sum((plot_mat[,1] - plot_mat[,2])^2/(plot_mat[,2] + plot_mat[,2]^2/x))/2
})

plot(r_seq, vec)
lines(c(-1e3, 1e3), rep(nrow(plot_mat), 2), col = "red", lwd = 2, lty = 2)
