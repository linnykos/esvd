rm(list=ls())
source("../experiment/Week24_simulation_generator.R")
load("../experiment/Week24_simulation_exponential.RData")

set.seed(10)
#simulation <- .data_generator(total = 200, col_drop = F)
#dat <- simulation$dat
# .plot_singlecell(dat)
length(which(dat == 0))/prod(dim(dat))
length(which(simulation$obs_mat2 == 0))/prod(dim(simulation$obs_mat2)) #percentage of true zeros

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #blue
             rgb(100,100,200,maxColorValue=255), #purple
             rgb(149,219,144,maxColorValue=255)) #green

png("../figure/experiment/25_latent.png", height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(simulation$cell_mat_org[,1], simulation$cell_mat_org[,2],
     xlim = range(c(simulation$cell_mat_org[,1], 0)),
     ylim = range(c(simulation$cell_mat_org[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

plot(simulation$gene_mat_org[,1], simulation$gene_mat_org[,2],
     xlim = range(c(simulation$gene_mat_org[,1], 0)),
     ylim = range(c(simulation$gene_mat_org[,2], 0)),
     col = col_vec[rep(1:simulation$g, each = simulation$d_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
graphics.off()

####

png("../figure/experiment/25_data.png", height = 1000, width = 2400, res = 300, units = "px")
par(mar = c(0.5, 0.5, 3, 0.5), mfrow = c(1,2))
.plot_singlecell(abs(simulation$gram_mat), main = "Inner product matrix")
lines(rep(0.5, 2), c(0,1), lwd = 5, lty = 2)
for(i in 1:3){
  lines(c(0, 1), rep(i/4, 2), lwd = 5, lty = 2)
}

.plot_singlecell(dat, main = "Observed matrix")
lines(rep(0.5, 2), c(0,1), lwd = 5, lty = 2)
for(i in 1:3){
  lines(c(0, 1), rep(i/4, 2), lwd = 5, lty = 2)
}
graphics.off()

###

png("../figure/experiment/25_true_zero.png", height = 1000, width = 1200, res = 300, units = "px")
par(mar = c(0.5, 0.5, 3, 0.5))
true_zero <- matrix(1, ncol = ncol(dat), nrow = nrow(dat))
true_zero[which(simulation$gram_mat < quantile(simulation$gram_mat, probs = 0.05))] <- 0
image(.rotate(true_zero), breaks = c(-0.5,0.5,1.5), col = c("blue3", "gray88"), asp = nrow(true_zero)/ncol(true_zero),
      axes = F, main = "Position of the true zeros")
lines(rep(0.5, 2), c(0,1), lwd = 5, lty = 2)
for(i in 1:3){
  lines(c(0, 1), rep(i/4, 2), lwd = 5, lty = 2)
}
graphics.off()

###

png("../figure/experiment/25_data_properties.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
zz <- apply(dat, 1, function(x){length(which(x != 0))/length(x)})
plot(sort(zz), ylab = "Percentage non-zero", main = "Sparsity per cell (row)",
     xlab = "Row rank", pch = 16, col = rgb(0,0,0,0.25))
zz <- apply(dat, 2, function(x){length(which(x != 0))/length(x)})
plot(sort(zz), ylab = "Percentage non-zero", main = "Sparsity per gene (colum)",
     xlab = "Column rank", pch = 16, col = rgb(0,0,0,0.25))
plot(sort(dat[dat != 0]), xlab = "Sorted rank (non-zeros)", ylab = "Value",
     main = "Non-zero quantiles", pch = 16, col = rgb(0,0,0,0.25))
graphics.off()

#####################

par(mfrow = c(1,3))
plot(res_nodropout$u_mat[,1], res_nodropout$u_mat[,2],
     xlim = range(c(res_nodropout$u_mat[,1], 0)),
     ylim = range(c(res_nodropout$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

plot(res_nodropout$v_mat[,1], res_nodropout$v_mat[,2],
     xlim = range(c(res_nodropout$v_mat[,1], 0)),
     ylim = range(c(res_nodropout$v_mat[,2], 0)),
     col = col_vec[rep(1:simulation$g, each = simulation$d_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

pred_mat <- res_nodropout$u_mat %*% t(res_nodropout$v_mat)
plot(as.numeric(simulation$gram_mat), as.numeric(pred_mat), asp = T)
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)

###

par(mfrow = c(1,2))
plot(res_withdropout$u_mat[,1], res_withdropout$u_mat[,2],
     xlim = range(c(res_withdropout$u_mat[,1], 0)),
     ylim = range(c(res_withdropout$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
pred_mat <- res_withdropout$u_mat %*% t(res_withdropout$v_mat)
plot(as.numeric(simulation$gram_mat), pred_mat, asp = T)
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)

par(mfrow = c(1,2))
plot(res_withimpute$u_mat[,1], res_withimpute$u_mat[,2],
     xlim = range(c(res_withimpute$u_mat[,1], 0)),
     ylim = range(c(res_withimpute$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
pred_mat <- res_withimpute$u_mat %*% t(res_withimpute$v_mat)
plot(as.numeric(simulation$gram_mat), pred_mat, asp = T)
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)

####################

zero_mat2 <- zero_mat
zero_mat2[is.na(zero_mat2)] <- 2
length(which(zero_mat2 == 0))/prod(dim(zero_mat2))

par(mar = rep(0.5,4))
image(.rotate(zero_mat2), breaks = c(-0.5,0.5,1.5,2.5), col = c("blue3", "gold", "firebrick1"), asp = nrow(zero_mat)/ncol(zero_mat),
      axes = F)

idx1 <- which(simulation$obs_mat < quantile(as.numeric(simulation$obs_mat[simulation$obs_mat > 0]), probs = 0.05))
idx2 <- which(zero_mat == 0)
length(intersect(idx1, idx2))/prod(dim(zero_mat))


