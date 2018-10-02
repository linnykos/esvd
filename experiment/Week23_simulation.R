rm(list=ls())
load("../experiment/Week23_simulation.RData")
source("../experiment/Week23_simulation_generator.R")

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #purple
             rgb(100,100,200,maxColorValue=255), #blue
             rgb(149,219,144,maxColorValue=255)) #green

###############

set.seed(10)
simulation2 <- .data_generator(total = 200)

library(slingshot)
slingshot_res <- slingshot::slingshot(data = simulation2$cell_mat_org,
                                      clusterLabels = rep(1:simulation2$h, each = simulation2$n_each),
                                      omega = 9)

png("../figure/experiment/23_latent.png", height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(simulation2$cell_mat_org[,1], simulation2$cell_mat_org[,2],
     xlim = range(c(simulation2$cell_mat_org[,1], 0)),
     ylim = range(c(simulation2$cell_mat_org[,2], 0)),
     col = col_vec[rep(1:simulation2$h, each = simulation2$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
lines(slingshot_res)

plot(simulation2$gene_mat_org[,1], simulation2$gene_mat_org[,2],
     xlim = range(c(simulation2$gene_mat_org[,1], 0)),
     ylim = range(c(simulation2$gene_mat_org[,2], 0)),
     col = col_vec[rep(1:simulation2$g, each = simulation2$d_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors")
lines(c(-1e6, 1e6), rep(0, 2), col = "red", lwd = 2, lty = 2)
lines( rep(0, 2), c(-1e6, 1e6), col = "red", lwd = 2, lty = 2)
graphics.off()

##################

png("../figure/experiment/23_data.png", height = 1000, width = 2400, res = 300, units = "px")
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

png("../figure/experiment/23_histograms.png", height = 2000, width = 2500, res = 300, units = "px")
par(mfrow = c(2,2), mar = c(4,3,3,1))
zz <- which.max(apply(dat, 1, function(x){length(which(x == 0))}))
.hist_augment(dat[zz,], main = paste0("Cell ", zz), xlab = "Value")
zz <- which.min(apply(dat, 1, function(x){length(which(x == 0))}))
.hist_augment(dat[zz,], main = paste0("Cell ", zz), xlab = "Value")
zz <- which.max(apply(dat, 2, function(x){length(which(x == 0))}))
.hist_augment(dat[,zz], main = paste0("Gene ", zz), xlab = "Value")
zz <- which.min(apply(dat, 2, function(x){length(which(x == 0))}))
.hist_augment(dat[,zz], main = paste0("Gene ", zz), xlab = "Value")
graphics.off()

###########################

# par(mfrow = c(1,3))
# zz <- apply(dat, 1, function(x){length(which(x != 0))/length(x)})
# plot(sort(zz))
# zz <- apply(dat, 2, function(x){length(which(x != 0))/length(x)})
# plot(sort(zz))
# plot(sort(dat[dat != 0]))


###########

# try the naive thing
dat_inv <- dat
dat_inv[which(dat > 0)] <- 1/dat[which(dat > 0)]
dat_inv[which(dat == 0)] <- max(dat_inv)

svd_res <- svd(dat_inv)
k <- 2
u_mat_naive <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
v_mat_naive <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
library(slingshot)
slingshot_res <- slingshot::slingshot(data = u_mat_naive, clusterLabels = rep(1:simulation$h, each = simulation$n_each))

plot(u_mat_naive[,1], u_mat_naive[,2],
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Estimated dim. 1", ylab = "Estimated dim. 2", main = "Cell estimated vectors")
lines(slingshot_res)

# let's see how well this naive method did
.evaluate_objective(dat, u_mat_naive, v_mat_naive)

##############

# compare to our implemented method
.evaluate_objective(dat, res_nodropout$u_mat, res_nodropout$v_mat)

plot(res_nodropout$u_mat[,1], res_nodropout$u_mat[,2],
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Estimated dim. 1", ylab = "Estimated dim. 2", main = "Cell estimated vectors")

