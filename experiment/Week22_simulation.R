rm(list=ls())
source("../experiment/Week22_simulation_generator.R")

col_vec <- c(rgb(205,40,54,maxColorValue=255),
             rgb(180,200,255,maxColorValue=255),
             rgb(100,100,200,maxColorValue=255),
             rgb(149,219,144,maxColorValue=255))

set.seed(10)
simulation <- .data_generator1(1, 2)
length(which(simulation$dat == 0))/prod(dim(simulation$dat))

library(slingshot)
slingshot_pop <- slingshot::slingshot(data = simulation$cell_mat,
                                      clusterLabels = rep(1:simulation$h,
                                                          each = simulation$n_each),
                                      omega = 8)

png("../figure/experiment/22_population.png", height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(simulation$cell_mat[,1], simulation$cell_mat[,2],
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], pch = 16,
     asp = T, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
lines(slingshot_pop)

plot(simulation$gene_mat[,1], simulation$gene_mat[,2],
     col = col_vec[rep(1:simulation$g, each = simulation$d_each)],
     asp = T, pch = 16,
     xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Gene latent vectors")
graphics.off()

png("../figure/experiment/22_observed.png", height = 1200, width = 1800, res = 300, units = "px")
par(mar = rep(0.5,4))
.plot_singlecell(simulation$dat)

for(i in 1:3){
  lines(c(0,1), rep(i/4, 2), lwd = 5, lty = 2)
}
lines(rep(1/2, 2), c(0,1), lwd = 5, lty = 2)
graphics.off()

png("../figure/experiment/22_observed_hist.png", height = 800, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,3,3,0.5))
yy <- apply(simulation$dat, 1, function(x){mean(x[x > 0])})
.hist_augment(simulation$dat[which.min(yy),], xlab = "Value", ylab = "Count",
              main = paste0("Cell ", which.min(yy)))
.hist_augment(simulation$dat[which.max(yy),], xlab = "Value", ylab = "Count",
              main = paste0("Cell ", which.max(yy)))
graphics.off()

###################

# now do the estimation
res_svd <- svd(simulation$dat)
u_mat <- res_svd$u[,1:simulation$k] %*% diag(sqrt(res_svd$d[1:simulation$k]))

slingshot_res <- slingshot::slingshot(data = u_mat,
                                      clusterLabels = rep(1:simulation$h, each = simulation$n_each),
                                      omega = 5)

png("../figure/experiment/22_estimated.png", height = 1200, width = 1000, res = 300, units = "px")
plot(u_mat[,1], u_mat[,2],
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Estimated dim. 1", ylab = "Estimated dim. 2", main = "Cell estimated vectors")
lines(slingshot_res)
graphics.off()

#########################

# hard setting
set.seed(10)
simulation <- .data_generator1(1, 2, distr_func = function(x){
  stats::rnorm(1, x, sd = abs(x))
})

png("../figure/experiment/22_observed_2.png", height = 1200, width = 1800, res = 300, units = "px")
par(mar = rep(0.5,4))
.plot_singlecell(simulation$dat)

for(i in 1:3){
  lines(c(0,1), rep(i/4, 2), lwd = 5, lty = 2)
}
lines(rep(1/2, 2), c(0,1), lwd = 5, lty = 2)
graphics.off()

length(which(simulation$dat == 0))/prod(dim(simulation$dat))

res_svd <- svd(simulation$dat)
u_mat <- res_svd$u[,1:simulation$k] %*% diag(sqrt(res_svd$d[1:simulation$k]))
slingshot_res <- slingshot::slingshot(data = u_mat,
                                      clusterLabels = rep(1:simulation$h, each = simulation$n_each))
png("../figure/experiment/22_estimated_2.png", height = 1200, width = 1000, res = 300, units = "px")
plot(u_mat[,1], u_mat[,2],
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Estimated dim. 1", ylab = "Estimated dim. 2", main = "Cell estimated vectors")
lines(slingshot_res)
graphics.off()

png("../figure/experiment/22_observed_hist_2.png", height = 800, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,3,3,0.5))
yy <- apply(simulation$dat, 1, function(x){mean(x[x > 0])})
.hist_augment(simulation$dat[which.min(yy),], xlab = "Value", ylab = "Count",
              main = paste0("Cell ", which.min(yy)))
.hist_augment(simulation$dat[which.max(yy),], xlab = "Value", ylab = "Count",
              main = paste0("Cell ", which.max(yy)))
graphics.off()

######

# different distribution
set.seed(10)
simulation <- .data_generator1(1,0, rotate_cells = T, distr_func = function(x){
  stats::rnorm(1, x, sd = abs(x))
})

png("../figure/experiment/22_population_3.png", height = 1200, width = 1000, res = 300, units = "px")
plot(simulation$cell_mat[,1], simulation$cell_mat[,2],
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], pch = 16,
     asp = T, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")
graphics.off()

png("../figure/experiment/22_observed_3.png", height = 1200, width = 1800, res = 300, units = "px")
par(mar = rep(0.5,4))
.plot_singlecell(simulation$dat)

for(i in 1:3){
  lines(c(0,1), rep(i/4, 2), lwd = 5, lty = 2)
}
lines(rep(1/2, 2), c(0,1), lwd = 5, lty = 2)
graphics.off()

length(which(simulation$dat == 0))/prod(dim(simulation$dat))

res_svd <- svd(simulation$dat)
u_mat <- res_svd$u[,1:simulation$k] %*% diag(sqrt(res_svd$d[1:simulation$k]))
slingshot_res <- slingshot::slingshot(data = u_mat,
                                      clusterLabels = rep(1:simulation$h, each = simulation$n_each))
png("../figure/experiment/22_estimated_3.png", height = 1200, width = 1000, res = 300, units = "px")
plot(u_mat[,1], u_mat[,2],
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Estimated dim. 1", ylab = "Estimated dim. 2", main = "Cell estimated vectors")
lines(slingshot_res)
graphics.off()

png("../figure/experiment/22_observed_hist_3.png", height = 800, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,3,3,0.5))
yy <- apply(simulation$dat, 1, function(x){mean(x[x > 0])})
.hist_augment(simulation$dat[which.min(yy),], xlab = "Value", ylab = "Count",
              main = paste0("Cell ", which.min(yy)))
.hist_augment(simulation$dat[which.max(yy),], xlab = "Value", ylab = "Count",
              main = paste0("Cell ", which.max(yy)))
graphics.off()

