rm(list=ls())
load("../results/factorization_esvd.RData")
# clean results
for(i in 1:length(res)){
  res[[i]] <- res[[i]][which(sapply(res[[i]], function(x){!all(is.na(x))}))]
}
sapply(res, function(x){
  length(x)
})

####################################

i <- 20
plot_mat <- lapply(1:length(res[[i]]), function(j){
  cbind(res[[i]][[j]]$missing_val, res[[i]][[j]]$expected_val)
})
if(length(plot_mat) > 1) plot_mat <- do.call(rbind, plot_mat) else plot_mat <- plot_mat[[1]]

set.seed(10)
new_sample <- sapply(plot_mat[,2], function(x){stats::rpois(1, lambda = x)})

plot(plot_mat[,2], new_sample, asp = T, pch = 16)
lines(c(-1e5,1e5), c(-1e5, 1e5), col = "red", lty = 2, lwd = 2)
seq_val <- seq(0, 4000, length.out = 500)
size <- as.numeric(paramMat[i,"true_r"])
y_bot <- sapply(seq_val, function(x){
  stats::qnbinom(0.1, size = size, prob = size/(size+x))
})
y_top <- sapply(seq_val, function(x){
  stats::qnbinom(0.9, size = size, prob = size/(size+x))
})
lines(seq_val, y_bot, col = "red", lty = 1, lwd = 2)
lines(seq_val, y_top, col = "red", lty = 1, lwd = 2)
pca_res <- stats::princomp(plot_mat)
lines(c(-1e10,1e10)*pca_res$loadings[1,2], c(-1e10,1e10)*pca_res$loadings[1,1],
      col = "blue", lwd = 2, lty = 2)

#########################

set.seed(10)
x = 4000
zz <- stats::rnbinom(10000, size = size, prob = 1-size/(size+x))
quantile(zz)

set.seed(10)
zz <- stats::rnbinom(10000, size = 50, prob = 0.8)
mean(zz)
50*(1-0.8)/0.8
50/(50+mean(zz))
