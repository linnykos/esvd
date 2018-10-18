rm(list=ls())
source("../experiment/Week27_simulation_generator.R")

# set.seed(10)
# simulation <- .data_generator_exponential(total = 100, col_drop = F)
# dat <- simulation$dat
# .plot_singlecell(dat)
length(which(dat == 0))/prod(dim(dat))
length(which(simulation$obs_mat2 == 0))/prod(dim(simulation$obs_mat2)) #percentage of true zeros

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #purple
             rgb(100,100,200,maxColorValue=255), #blue
             rgb(149,219,144,maxColorValue=255)) #green

png("../figure/experiment/25_latent.png", height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(simulation$cell_mat[,1], simulation$cell_mat[,2],
     xlim = range(c(simulation$cell_mat[,1], 0)),
     ylim = range(c(simulation$cell_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, 1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, 1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)

plot(simulation$gene_mat[,1], simulation$gene_mat[,2],
     xlim = range(c(simulation$gene_mat[,1], 0)),
     ylim = range(c(simulation$gene_mat[,2], 0)),
     col = col_vec[rep(1:simulation$g, each = simulation$d_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, -1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, -1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)

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

####

png("../figure/experiment/25_true_scatter.png", height = 1000, width = 1000, res = 300, units = "px")
plot(-1/simulation$gram_mat, simulation$obs_mat, pch = 16, col = rgb(0,0,0,0.1),
     xlim = c(0,1.2), asp = T,
     ylab = "True observed value", xlab = "True expected value",
     main = "Simulated data")
lines(c(-1e6, 1e6), c(-1e6, 1e6), col = "red", lty = 2, lwd = 2)
lines(rep(0,2), c(-1e6, 1e6), col = "red", lty = 2, lwd = 1)
lines(c(-1e6, 1e6), rep(0,2), col = "red", lty = 2, lwd = 1)

# plot confidence bands
zz <- seq(0, 5, length.out = 500)
upper <- sapply(zz, function(x){qexp(p = 0.9, rate = 1/x)})
lower <- sapply(zz, function(x){qexp(p = 0.1, rate = 1/x)})

lines(zz, upper, col = "red", lty = "dotted", lwd = 1)
lines(zz, lower, col = "red", lty = "dotted", lwd = 1)

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

# IDEAL FIT
png("../figure/experiment/25_fit_ideal.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(res_ideal$u_mat[,1], res_ideal$u_mat[,2],
     xlim = range(c(res_ideal$u_mat[,1], 0)),
     ylim = range(c(res_ideal$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(Ideal fit)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, 1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, 1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)

plot(res_ideal$v_mat[,1], res_ideal$v_mat[,2],
     xlim = range(c(res_ideal$v_mat[,1], 0)),
     ylim = range(c(res_ideal$v_mat[,2], 0)),
     col = col_vec[rep(1:simulation$g, each = simulation$d_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors\n(Ideal fit)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, -1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, -1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)

pred_mat <- res_ideal$u_mat %*% t(res_ideal$v_mat)
plot(as.numeric(simulation$gram_mat), as.numeric(pred_mat),
     pch = 16, col = rgb(0,0,0,0.1), asp = T,
     xlab = "True inner product value",
     ylab = "Predicted inner product value",
     main = "Gram matrix\n(Ideal fit)")
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)

err <- mean((simulation$gram_mat - pred_mat)^2)
coords <- par("usr")
text(x = .05*coords[1]+0.95*coords[2], y = .9*coords[3]+0.1*coords[4],
     labels = paste0("Error: ", round(err, 2)),
     col = "red", pos = 2)
graphics.off()

###

#NO DROPOUT
png("../figure/experiment/25_fit_nodropout.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(res_nodropout$u_mat[,1], res_nodropout$u_mat[,2],
     xlim = range(c(res_nodropout$u_mat[,1], 0)),
     ylim = range(c(res_nodropout$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(No dropout)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, 1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, 1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)

plot(res_nodropout$v_mat[,1], res_nodropout$v_mat[,2],
     xlim = range(c(res_nodropout$v_mat[,1], 0)),
     ylim = range(c(res_nodropout$v_mat[,2], 0)),
     col = col_vec[rep(1:simulation$g, each = simulation$d_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors\n(No dropout)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, -1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, -1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)

pred_mat <- res_nodropout$u_mat %*% t(res_nodropout$v_mat)
plot(as.numeric(simulation$gram_mat), as.numeric(pred_mat),
     pch = 16, col = rgb(0,0,0,0.1), asp = T,
     xlab = "True inner product value",
     ylab = "Predicted inner product value",
     main = "Gram matrix\n(No dropout)")
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)

err <- mean((simulation$gram_mat - pred_mat)^2)
coords <- par("usr")
text(x = .05*coords[1]+0.95*coords[2], y = .9*coords[3]+0.1*coords[4],
     labels = paste0("Error: ", round(err, 2)),
     col = "red", pos = 2)
graphics.off()

#NO DROPOUT BUT CHEAT
png("../figure/experiment/25_fit_nodropout_cheat.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(res_nodropout_cheat$u_mat[,1], res_nodropout_cheat$u_mat[,2],
     xlim = range(c(res_nodropout_cheat$u_mat[,1], 0)),
     ylim = range(c(res_nodropout_cheat$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(No dropout, with cheating)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, 1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, 1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)

plot(res_nodropout_cheat$v_mat[,1], res_nodropout_cheat$v_mat[,2],
     xlim = range(c(res_nodropout_cheat$v_mat[,1], 0)),
     ylim = range(c(res_nodropout_cheat$v_mat[,2], 0)),
     col = col_vec[rep(1:simulation$g, each = simulation$d_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors\n(No dropout, with cheating)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, -1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, -1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)

pred_mat <- res_nodropout_cheat$u_mat %*% t(res_nodropout_cheat$v_mat)
plot(as.numeric(simulation$gram_mat), as.numeric(pred_mat),
     pch = 16, col = rgb(0,0,0,0.1), asp = T,
     xlab = "True inner product value",
     ylab = "Predicted inner product value",
     main = "Gram matrix\n(No dropout, with cheating)")
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)

err <- mean((simulation$gram_mat - pred_mat)^2)
coords <- par("usr")
text(x = .05*coords[1]+0.95*coords[2], y = .9*coords[3]+0.1*coords[4],
     labels = paste0("Error: ", round(err, 2)),
     col = "red", pos = 2)
graphics.off()

res_nodropout$obj_vec[length(res_nodropout$obj_vec)]
res_nodropout_cheat$obj_vec[length(res_nodropout_cheat$obj_vec)]

#########################

# INVESTIGATION INTO IMPUTATION
# true zeros

png("../figure/experiment/25_true_zero2.png", height = 1000, width = 2400, res = 300, units = "px")
par(mar = c(0.5, 0.5, 3, 1), mfrow = c(1,2))
true_zero <- matrix(1, ncol = ncol(dat), nrow = nrow(dat))
true_zero[which(simulation$gram_mat < quantile(simulation$gram_mat, probs = 0.05))] <- 0
image(.rotate(true_zero), breaks = c(-0.5,0.5,1.5), col = c("blue3", "gray88"), asp = nrow(true_zero)/ncol(true_zero),
      axes = F, main = "Position of the true zeros")
lines(rep(0.5, 2), c(0,1), lwd = 5, lty = 2)
for(i in 1:3){
  lines(c(0, 1), rep(i/4, 2), lwd = 5, lty = 2)
}

est_zero <- zero_mat
est_zero[is.na(est_zero)] <- 2
image(.rotate(est_zero), breaks = c(-0.5,0.5,1.5,2.5), col = c("blue3", "gray88", "goldenrod1"), asp = nrow(true_zero)/ncol(true_zero),
      axes = F, main = "Position of the estimated zeros")
lines(rep(0.5, 2), c(0,1), lwd = 5, lty = 2)
for(i in 1:3){
  lines(c(0, 1), rep(i/4, 2), lwd = 5, lty = 2)
}
graphics.off()

length(which(zero_mat == 0))/prod(dim(zero_mat))

# histograms on to-be imputed values
png("../figure/experiment/25_hist_values.png", height = 800, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
max_val <- 0.15
idx <- which(zero_mat == 0)
.hist_augment(dat[idx], max_val = max_val, xlab = "Value", main = "Histogram of true zeros")

idx <- which(is.na(zero_mat))
.hist_augment(dat[idx], max_val = max_val, xlab = "Value", main = "Histogram of to-be imputed values")

idx <- which(zero_mat == 1)
.hist_augment(dat[idx], max_val = max_val, xlab = "Value", main = "Histogram of untounched values")
graphics.off()

# now for the values we imputed to
png("../figure/experiment/25_imputed_scatter.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(5, 5, 0.5, 1))
idx <- which(is.na(zero_mat))
# plot(dat_impute[idx], dat[idx], asp = T)
plot(dat_impute[idx], simulation$obs_mat[idx], asp = T, pch = 16, col = rgb(0,0,0,0.1),
     xlab = "Imputed value", ylab = "True observation value",
     ylim = c(0, 0.3))
lines(c(0, 1e6), c(0, 1e6), col = "red", lty = 2, lwd = 2)
lines(rep(0,2), c(-1e6, 1e6), col = "red", lty = 2, lwd = 1)
lines(c(-1e6, 1e6), rep(0,2), col = "red", lty = 2, lwd = 1)
err <- sum((dat_impute[idx] - simulation$obs_mat[idx])^2)
coords <- par("usr")
text(x = .05*coords[1]+0.95*coords[2], y = .1*coords[3]+0.9*coords[4],
     labels = paste0("Error: ", round(err, 2)),
     col = "red", pos = 2)

plot(dat_impute[idx], -1/simulation$gram_mat[idx], asp = T, pch = 16, col = rgb(0,0,0,0.1),
     xlab = "Imputed value", ylab = "True expected value",
     ylim = c(0, 0.3))
lines(c(-1e6, 1e6), c(-1e6, 1e6), col = "red", lty = 2, lwd = 2)
lines(rep(0,2), c(-1e6, 1e6), col = "red", lty = 2, lwd = 1)
lines(c(-1e6, 1e6), rep(0,2), col = "red", lty = 2, lwd = 1)
err <- sum((dat_impute[idx] - -1/simulation$gram_mat[idx])^2)
coords <- par("usr")
text(x = .05*coords[1]+0.95*coords[2], y = .1*coords[3]+0.9*coords[4],
     labels = paste0("Error: ", round(err, 2)),
     col = "red", pos = 2)
graphics.off()


############################
# now for our usual histograms

png("../figure/experiment/25_fit_withdropout.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(res_withdropout$u_mat[,1], -res_withdropout$u_mat[,2],
     xlim = range(c(res_withdropout$u_mat[,1], 0)),
     ylim = range(c(-res_withdropout$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(With dropout)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, 1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, 1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)

plot(res_withdropout$v_mat[,1], -res_withdropout$v_mat[,2],
     xlim = range(c(res_withdropout$v_mat[,1], 0)),
     ylim = range(c(-res_withdropout$v_mat[,2], 0)),
     col = col_vec[rep(1:simulation$g, each = simulation$d_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors\n(With dropout)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, -1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, -1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)

pred_mat <- res_withdropout$u_mat %*% t(res_withdropout$v_mat)
plot(as.numeric(simulation$gram_mat), as.numeric(pred_mat), asp = T,
     pch = 16, col = rgb(0,0,0,0.1),
     xlab = "True inner product value",
     ylab = "Predicted inner product value",
     main = "Gram matrix\n(With dropout)")
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)

err <- mean((simulation$gram_mat - pred_mat)^2)
coords <- par("usr")
text(x = .05*coords[1]+0.95*coords[2], y = .9*coords[3]+0.1*coords[4],
     labels = paste0("Error: ", round(err, 2)),
     col = "red", pos = 2)
graphics.off()

png("../figure/experiment/25_fit_withdropout_cheat.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(res_withdropout_cheat$u_mat[,1], res_withdropout_cheat$u_mat[,2],
     xlim = range(c(res_withdropout_cheat$u_mat[,1], 0)),
     ylim = range(c(res_withdropout_cheat$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(With dropout, with cheating)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, 1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, 1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)

plot(res_withdropout_cheat$v_mat[,1], res_withdropout_cheat$v_mat[,2],
     xlim = range(c(res_withdropout_cheat$v_mat[,1], 0)),
     ylim = range(c(res_withdropout_cheat$v_mat[,2], 0)),
     col = col_vec[rep(1:simulation$g, each = simulation$d_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors\n(With dropout, with cheating)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, -1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, -1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)

pred_mat <- res_withdropout_cheat$u_mat %*% t(res_withdropout_cheat$v_mat)
plot(as.numeric(simulation$gram_mat), as.numeric(pred_mat), asp = T,
     pch = 16, col = rgb(0,0,0,0.1),
     xlab = "True inner product value",
     ylab = "Predicted inner product value",
     main = "Gram matrix\n(With dropout, with cheating)")
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)

err <- mean((simulation$gram_mat - pred_mat)^2)
coords <- par("usr")
text(x = .05*coords[1]+0.95*coords[2], y = .9*coords[3]+0.1*coords[4],
     labels = paste0("Error: ", round(err, 2)),
     col = "red", pos = 2)
graphics.off()

###############################
png("../figure/experiment/25_fit_withimpute.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(res_withimpute$u_mat[,1], res_withimpute$u_mat[,2],
     xlim = range(c(res_withimpute$u_mat[,1], 0)),
     ylim = range(c(res_withimpute$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(With impute)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, 1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, 1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)

plot(res_withimpute$v_mat[,1], res_withimpute$v_mat[,2],
     xlim = range(c(res_withimpute$v_mat[,1], 0)),
     ylim = range(c(res_withimpute$v_mat[,2], 0)),
     col = col_vec[rep(1:simulation$g, each = simulation$d_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors\n(With impute)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, -1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, -1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)

pred_mat <- res_withimpute$u_mat %*% t(res_withimpute$v_mat)
plot(as.numeric(simulation$gram_mat), as.numeric(pred_mat), asp = T,
     pch = 16, col = rgb(0,0,0,0.1),
     xlab = "True inner product value",
     ylab = "Predicted inner product value",
     main = "Gram matrix\n(With impute)")
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)

err <- mean((simulation$gram_mat - pred_mat)^2)
coords <- par("usr")
text(x = .05*coords[1]+0.95*coords[2], y = .9*coords[3]+0.1*coords[4],
     labels = paste0("Error: ", round(err, 2)),
     col = "red", pos = 2)
graphics.off()

png("../figure/experiment/25_fit_withimpute_cheat.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(res_withimpute_cheat$u_mat[,1], res_withimpute_cheat$u_mat[,2],
     xlim = range(c(res_withimpute_cheat$u_mat[,1], 0)),
     ylim = range(c(res_withimpute_cheat$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(With impute, with cheating)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, 1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, 1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)

plot(res_withimpute_cheat$v_mat[,1], res_withimpute_cheat$v_mat[,2],
     xlim = range(c(res_withimpute_cheat$v_mat[,1], 0)),
     ylim = range(c(res_withimpute_cheat$v_mat[,2], 0)),
     col = col_vec[rep(1:simulation$g, each = simulation$d_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors\n(With impute, with cheating)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, -1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, -1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)

pred_mat <- res_withimpute_cheat$u_mat %*% t(res_withimpute_cheat$v_mat)
plot(as.numeric(simulation$gram_mat), as.numeric(pred_mat), asp = T,
     pch = 16, col = rgb(0,0,0,0.1),
     xlab = "True inner product value",
     ylab = "Predicted inner product value",
     main = "Gram matrix\n(With impute, with cheating)")
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)

err <- mean((simulation$gram_mat - pred_mat)^2)
coords <- par("usr")
text(x = .05*coords[1]+0.95*coords[2], y = .9*coords[3]+0.1*coords[4],
     labels = paste0("Error: ", round(err, 2)),
     col = "red", pos = 2)
graphics.off()

####################

# do the "naive" analyses
svd_res <- svd(dat)
k <- 6
u_mat_naive <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
v_mat_naive <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))

set.seed(10)
svd_impute <- singlecell:::.initialization(dat_impute, family = "exponential",
                                           max_val = max_val)

png("../figure/experiment/25_naive.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(u_mat_naive[,1], u_mat_naive[,2],
     xlim = range(c(u_mat_naive[,1], 0)),
     ylim = range(c(u_mat_naive[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(Naive fit, mean)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

plot(init$u_mat[,1], init$u_mat[,2],
     xlim = range(c(init$u_mat[,1], 0)),
     ylim = range(c(init$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(Naive fit, -1/mean)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

plot(svd_impute$u_mat[,1], svd_impute$u_mat[,2],
     xlim = range(c(svd_impute$u_mat[,1], 0)),
     ylim = range(c(svd_impute$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(Naive fit, -1/(imputed mean))")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

graphics.off()

