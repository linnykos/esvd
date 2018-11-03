rm(list=ls())
source("../experiment/Week27_simulation_generator.R")
load("../experiment/Week27_simulation.RData")

# ###################
# # fancy way to rearrange the rows
# row_sparsity <- apply(dat, 1, function(x){sum(x != 0)})
#
# # rearrange within each 50
# row_idx_vec <- rep(NA, nrow(dat))
# width <- 4
# for(i in 1:4){
#   tmp_vec <- ((i-1)*50+1):(i*50)
#   ord <- order(row_sparsity[tmp_vec], decreasing = T)
#
#   tmp_mat <- cbind(ord[seq(1, length(ord), by = 2)], ord[seq(2, length(ord), by = 2)])
#   ord <- as.numeric(tmp_mat)
#
#   for(j in (width+1):(length(tmp_vec)-width-1)){
#     tmp_idx <- (j-width):(j+width)
#     ord[tmp_idx] <- ord[sample(tmp_idx)]
#   }
#
#   row_idx_vec[tmp_vec] <- ord + (i-1)*50
# }
#
# row_sparsity <- apply(dat, 2, function(x){sum(x != 0)})
# col_idx_vec <- rep(NA, ncol(dat))
# width <- 4
# for(i in 1:2){
#   tmp_vec <- ((i-1)*120+1):(i*120)
#   ord <- order(row_sparsity[tmp_vec], decreasing = T)
#
#   tmp_mat <- cbind(ord[seq(1, length(ord), by = 2)], ord[seq(2, length(ord), by = 2)])
#   ord <- as.numeric(tmp_mat)
#
#   for(j in (width+1):(length(tmp_vec)-width-1)){
#     tmp_idx <- (j-width):(j+width)
#     ord[tmp_idx] <- ord[sample(tmp_idx)]
#   }
#
#   col_idx_vec[tmp_vec] <- ord + (i-1)*120
# }
#
# dat <- dat[row_idx_vec,col_idx_vec]

####################

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #purple
             rgb(100,100,200,maxColorValue=255), #blue
             rgb(149,219,144,maxColorValue=255)) #green

cluster_labels <- rep(1:simulation$h, each = simulation$n_each)

res <- .get_lineages(simulation$cell_mat, cluster_labels)
lineages <- res$lineages
lineages[[1]] <- c("1","3","2")
lineages[[2]] <- c("1","3","4")
cluster_mat <- res$cluster_mat
curves <- .get_curves_tmp(simulation$cell_mat, cluster_mat, lineages, reassign = F)


png("../figure/experiment/27_latent.png", height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(simulation$cell_mat[,1], simulation$cell_mat[,2],
     xlim = range(c(simulation$cell_mat[,1], 0)),
     ylim = range(c(simulation$cell_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")

for(i in 1:length(curves)){
  ord <- curves[[i]]$ord
  lines(curves[[i]]$s[ord,1], curves[[i]]$s[ord,2], lwd = 2)
}

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


png("../figure/experiment/27_inner_product.png", height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(0.5, 0.5, 3, 0.5))

tmp <- simulation$gram_mat - min(simulation$gram_mat)

.plot_singlecell(tmp, main = "Inner product matrix")
lines(rep(0.5, 2), c(0,1), lwd = 5, lty = 2)
for(i in 1:3){
  lines(c(0, 1), rep(i/4, 2), lwd = 5, lty = 2)
}
graphics.off()

png("../figure/experiment/27_data_uncorrupted.png", height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(0.5, 0.5, 3, 0.5))

.plot_singlecell(abs(simulation$obs_mat), main = "Fully-observed matrix", luminosity = F)
lines(rep(0.5, 2), c(0,1), lwd = 5, lty = 2)
for(i in 1:3){
  lines(c(0, 1), rep(i/4, 2), lwd = 5, lty = 2)
}
graphics.off()

png("../figure/experiment/27_data_uncorrupted_withzero.png", height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(0.5, 0.5, 3, 0.5))

dat_tmp <- simulation$obs_mat
dat_tmp[which(simulation$gram_mat < quantile(simulation$gram_mat, probs = 0.05))] <- 0
.plot_singlecell(dat_tmp, main = "Fully-observed matrix", luminosity = F)
lines(rep(0.5, 2), c(0,1), lwd = 5, lty = 2)
for(i in 1:3){
  lines(c(0, 1), rep(i/4, 2), lwd = 5, lty = 2)
}
graphics.off()




png("../figure/experiment/27_data.png", height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(0.5, 0.5, 3, 0.5))

.plot_singlecell(dat, main = "Observed matrix")
lines(rep(0.5, 2), c(0,1), lwd = 5, lty = 2)
for(i in 1:3){
  lines(c(0, 1), rep(i/4, 2), lwd = 5, lty = 2)
}
graphics.off()

###

png("../figure/experiment/27_true_zero.png", height = 1000, width = 1200, res = 300, units = "px")
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

#####################

# IDEAL FIT
curves <- .get_curves_tmp(res_ideal$u_mat, cluster_mat, lineages, reassign = F,
                          b = 3)

png("../figure/experiment/27_fit_ideal.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(res_ideal$u_mat[,1], res_ideal$u_mat[,2],
     xlim = range(c(res_ideal$u_mat[,1], 0)),
     ylim = range(c(res_ideal$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(Ideal fit)")

for(i in 1:length(curves)){
  ord <- curves[[i]]$ord
  lines(curves[[i]]$s[ord,1], curves[[i]]$s[ord,2], lwd = 2)
}

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

####################

# INVESTIGATION INTO IMPUTATION
# true zeros

png("../figure/experiment/27_true_zero2.png", height = 1000, width = 2400, res = 300, units = "px")
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

png("../figure/experiment/27_true_zero3.png", height = 1000, width = 2400, res = 300, units = "px")
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
image(.rotate(est_zero), breaks = c(-0.5,0.5,1.5,2.5), col = c("blue3", "gray88", "gray88"), asp = nrow(true_zero)/ncol(true_zero),
      axes = F, main = "Position of the estimated zeros")
lines(rep(0.5, 2), c(0,1), lwd = 5, lty = 2)
for(i in 1:3){
  lines(c(0, 1), rep(i/4, 2), lwd = 5, lty = 2)
}
graphics.off()

png("../figure/experiment/27_true_zero_imputed.png", height = 1000, width = 2400, res = 300, units = "px")
par(mar = c(0.5, 0.5, 3, 1), mfrow = c(1,2))
.plot_singlecell(dat_impute2, main = "Data after imputation", luminosity = F)
lines(rep(0.5, 2), c(0,1), lwd = 5, lty = 2)
for(i in 1:3){
  lines(c(0, 1), rep(i/4, 2), lwd = 5, lty = 2)
}

est_zero <- zero_mat
est_zero[is.na(est_zero)] <- 2
image(.rotate(est_zero), breaks = c(-0.5,0.5,1.5,2.5), col = c("blue3", "gray88", "goldenrod1"),
      asp = nrow(zero_mat)/ncol(zero_mat),
      axes = F, main = "Position of the estimated zeros")
lines(rep(0.5, 2), c(0,1), lwd = 5, lty = 2)
for(i in 1:3){
  lines(c(0, 1), rep(i/4, 2), lwd = 5, lty = 2)
}
graphics.off()


# histograms on to-be imputed values
png("../figure/experiment/27_hist_values.png", height = 800, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
max_val <- 0.15
idx <- which(zero_mat == 0)
.hist_augment(dat[idx], max_val = max_val, xlab = "Value", main = "Histogram of true zeros")

idx <- which(is.na(zero_mat))
.hist_augment(dat[idx], max_val = max_val, xlab = "Value", main = "Histogram of to-be imputed values")

idx <- which(zero_mat == 1)
.hist_augment(dat[idx], max_val = max_val, xlab = "Value", main = "Histogram of untounched values")
graphics.off()

############################
png("../figure/experiment/27_fit_nodropout.png", height = 1000, width = 2400, res = 300, units = "px")
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
plot(as.numeric(simulation$gram_mat), as.numeric(pred_mat), asp = T,
     pch = 16, col = rgb(0,0,0,0.1),
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

############################
curves <- .get_curves_tmp(res_withimpute2$u_mat, cluster_mat, lineages, reassign = F,
                          b = 4)

png("../figure/experiment/27_fit_withimpute2.png", height = 1000, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(res_withimpute2$u_mat[,1], res_withimpute2$u_mat[,2],
     xlim = range(c(res_withimpute2$u_mat[,1], 0)),
     ylim = range(c(res_withimpute2$u_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(With impute)")

for(i in 1:length(curves)){
  ord <- curves[[i]]$ord
  lines(curves[[i]]$s[ord,1], curves[[i]]$s[ord,2], lwd = 2)
}

lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, 1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, 1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)

plot(res_withimpute2$v_mat[,1], res_withimpute2$v_mat[,2],
     xlim = range(c(res_withimpute2$v_mat[,1], 0)),
     ylim = range(c(res_withimpute2$v_mat[,2], 0)),
     col = col_vec[rep(1:simulation$g, each = simulation$d_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors\n(With impute)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(c(0, -1e6), c(0, -1e6), col = "red", lwd = 1, lty = 2)
lines(c(0, -1e6), c(0, 1e6), col = "red", lwd = 1, lty = 2)

pred_mat <- res_withimpute2$u_mat %*% t(res_withimpute2$v_mat)
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


##########################

png("../figure/experiment/27_imputed_scatter.png", height = 1000, width = 2400, res = 300, units = "px")
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

####################

# do the "naive" analyses
svd_res <- svd(dat_impute)
k <- 6
u_mat_naive <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
v_mat_naive <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))

pred_mat <- init_impute$u_mat %*% t(init_impute$v_mat)
svd_res <- svd(pred_mat)
u_mat2 <- svd_res$u[,1:2] %*% diag(sqrt(svd_res$d[1:2]))

png("../figure/experiment/27_naive.png", height = 1000, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(u_mat_naive[,1], u_mat_naive[,2],
     xlim = range(c(u_mat_naive[,1], 0)),
     ylim = range(c(u_mat_naive[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(Naive fit, imputed mean)")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

plot(u_mat2[,1], u_mat2[,2],
     xlim = range(c(u_mat2[,1], 0)),
     ylim = range(c(u_mat2[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors\n(Naive fit, -1/(imputed mean))")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)

graphics.off()

