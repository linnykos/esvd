rm(list=ls())
source("../experiment/Week24_simulation_generator.R")
load("../experiment/Week24_simulation.RData")

set.seed(10)
#simulation <- .data_generator(total = 200, col_drop = F)
#dat <- simulation$dat
.plot_singlecell(dat)
length(which(dat == 0))/prod(dim(dat))
length(which(simulation$obs_mat2 == 0))/prod(dim(simulation$obs_mat2)) #percentage of true zeros

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #blue
             rgb(100,100,200,maxColorValue=255), #purple
             rgb(149,219,144,maxColorValue=255)) #green

plot(simulation$cell_mat[,1], simulation$cell_mat[,2],
     xlim = range(c(simulation$cell_mat[,1], 0)),
     ylim = range(c(simulation$cell_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

##################

par(mfrow = c(1,2))
plot(res_nodropout$u_mat[,1], res_nodropout$u_mat[,2],
     xlim = range(c(res_nodropout$u_mat[,1], 0)),
     ylim = range(c(res_nodropout$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
pred_mat <- res_nodropout$u_mat %*% t(res_nodropout$v_mat)
plot(as.numeric(simulation$gram_mat), pred_mat, asp = T)
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)

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
