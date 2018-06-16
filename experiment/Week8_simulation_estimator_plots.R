rm(list=ls())

# set up parameters
adj <- matrix(c(0, 0.1, -.05,
                0.15,-.1,-.1,
                -.25,-.15,0.05), 3, 3, byrow = T)
tmp <- svd(adj)
u_center <- t(tmp$u %*% diag(sqrt(tmp$d)))
v_center <- t(tmp$v %*% diag(sqrt(tmp$d)))

t(u_center) %*% v_center

u_num <- c(20, 10, 80)
u_label <- unlist(lapply(1:3, function(x){rep(x, u_num[x])}))
v_num <- c(30, 40, 100)
v_label <- unlist(lapply(1:3, function(x){rep(x, v_num[x])}))

# generate matrices
set.seed(10)
u_sig <- 0.03
u_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = u_num[x], mu = u_center[,x], Sigma = u_sig*diag(3))
}))
v_sig <- 0.03
v_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = v_num[x], mu = v_center[,x], Sigma = v_sig*diag(3))
}))

mean_dat <- u_dat %*% t(v_dat)
n <- nrow(mean_dat)
d <- ncol(mean_dat)
dat <- matrix(0, n, d)
setting <- matrix(0, n, d)
nodropout_dat <- matrix(0, n, d)

dropout_create <- function(a, b){
  function(x){ 1/(1+exp(-(a+b*x))) }
}
dropout_coef <- c(-1, 5)
dropout_function <- dropout_create(dropout_coef[1], dropout_coef[2])

set.seed(10)
for(i in 1:n){
  for(j in 1:d){
    if(mean_dat[i,j] <= 0) {
      dat[i,j] <- 0
      setting[i,j] <- 0
      nodropout_dat[i,j] <- 0
    } else {
      val <- rexp(1, rate = 1/mean_dat[i,j])
      bool <- rbinom(1, 1, prob = dropout_function(val))
      setting[i,j] <- bool + 1
      dat[i,j] <- bool*val
      nodropout_dat[i,j] <- mean_dat[i,j]
    }
  }
}

#####################

# our estimator
res <- estimate_latent(dat, k = 3, dropout_function, threshold = 0.3, verbose = T,
                       initialization = "SVD", alpha = 1000)

# naive estimator
res_svd <- svd(dat)
k <- 3
u_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

# truth
tmp <- u_dat %*% t(v_dat)
res_svd <- svd(tmp)
u_dat_rescaled <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_dat_rescaled <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

#######################

# plot of eigenspaces

png(paste0("../figure/experiment/8_simulation_cell_eigenspace.png"), height = 900, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(NA, xlim = range(u_dat_rescaled[,1]), ylim = range(u_dat_rescaled[,2]), asp = T,
     xlab = "U[,1]", ylab = "U[,2]", main = "True cell eigenspace")
lines(c(-10,10), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-10,10), col = "red", lwd = 2, lty = 2)
points(u_dat_rescaled[,1], u_dat_rescaled[,2], pch = 16, col = rgb(0,0,0,0.5))

plot(NA, xlim = range(res$u_mat[,1]), ylim = range(res$u_mat[,2]), asp = T,
     xlab = "U[,1]", ylab = "U[,2]", main = "Estimated cell eigenspace")
lines(c(-10,10), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-10,10), col = "red", lwd = 2, lty = 2)
points(res$u_mat[,1], res$u_mat[,2], pch = 16, col = rgb(0,0,0,0.5))

plot(NA, xlim = range(u_mat[,1]), ylim = range(u_mat[,2]), asp = T,
     xlab = "U[,1]", ylab = "U[,2]", main = "Naive estimated cell eigenspace")
lines(c(-10,10), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-10,10), col = "red", lwd = 2, lty = 2)
points(u_mat[,1], u_mat[,2], pch = 16, col = rgb(0,0,0,0.5))
graphics.off()

png(paste0("../figure/experiment/8_simulation_gene_eigenspace.png"), height = 900, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(NA, xlim = range(u_dat_rescaled[,1]), ylim = range(u_dat_rescaled[,2]), asp = T,
     xlab = "V[,1]", ylab = "V[,2]", main = "True gene eigenspace")
lines(c(-10,10), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-10,10), col = "red", lwd = 2, lty = 2)
points(u_dat_rescaled[,1], u_dat_rescaled[,2], pch = 16, col = rgb(0,0,0,0.5))

plot(NA, xlim = range(res$u_mat[,1]), ylim = range(res$u_mat[,2]), asp = T,
     xlab = "V[,1]", ylab = "V[,2]", main = "Estimated gene eigenspace")
lines(c(-10,10), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-10,10), col = "red", lwd = 2, lty = 2)
points(res$u_mat[,1], res$u_mat[,2], pch = 16, col = rgb(0,0,0,0.5))

plot(NA, xlim = range(u_mat[,1]), ylim = range(u_mat[,2]), asp = T,
     xlab = "V[,1]", ylab = "V[,2]", main = "Naive estimated gene eigenspace")
lines(c(-10,10), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-10,10), col = "red", lwd = 2, lty = 2)
points(u_mat[,1], u_mat[,2], pch = 16, col = rgb(0,0,0,0.5))
graphics.off()

############

png(paste0("../figure/experiment/8_simulation_cell_eigenspace_comparison.png"), height = 1100, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(NA, xlim = range(u_dat_rescaled[,1]), ylim = range(res$u_mat[,1]), asp = T,
     xlab = "U[,1] of True", ylab = "U[,1] of estimate", main = "Cell eigenspace comparison\nagainst our estimator")
lines(c(-10,10), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-10,10), col = "red", lwd = 2, lty = 2)
lines(c(-10,10), c(-10,10), col = "red", lwd = 2)
points(u_dat_rescaled[,1], res$u_mat[,1], pch = 16, col = rgb(0,0,0,0.5))

plot(NA, xlim = range(u_dat_rescaled[,1]), ylim = range(-u_mat[,1]), asp = T,
     xlab = "U[,1] of True", ylab = "U[,1] of estimate", main = "Cell eigenspace comparison\nagainst naive estimator")
lines(c(-10,10), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-10,10), col = "red", lwd = 2, lty = 2)
lines(c(-10,10), c(-10,10), col = "red", lwd = 2)
points(u_dat_rescaled[,1], -u_mat[,1], pch = 16, col = rgb(0,0,0,0.5))
graphics.off()

png(paste0("../figure/experiment/8_simulation_gene_eigenspace_comparison.png"), height = 1100, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(NA, xlim = range(v_dat_rescaled[,1]), ylim = range(res$v_mat[,1]), asp = T,
     xlab = "V[,1] of True", ylab = "V[,1] of estimate", main = "Gene eigenspace comparison\nagainst our estimator")
lines(c(-10,10), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-10,10), col = "red", lwd = 2, lty = 2)
lines(c(-10,10), c(-10,10), col = "red", lwd = 2)
points(v_dat_rescaled[,1], res$v_mat[,1], pch = 16, col = rgb(0,0,0,0.5))

plot(NA, xlim = range(v_dat_rescaled[,1]), ylim = range(-v_mat[,1]), asp = T,
     xlab = "V[,1] of True", ylab = "V[,1] of estimate", main = "Gene eigenspace comparison\nagainst naive estimator")
lines(c(-10,10), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-10,10), col = "red", lwd = 2, lty = 2)
lines(c(-10,10), c(-10,10), col = "red", lwd = 2)
points(v_dat_rescaled[,1], -v_mat[,1], pch = 16, col = rgb(0,0,0,0.5))
graphics.off()

###################### #predicted values plot
res_pred <- res$u_mat %*% t(res$v_mat)
res_pred_naive <- u_mat %*% t(v_mat)

png(paste0("../figure/experiment/8_simulation_prediction.png"), height = 1100, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(NA, xlim = range(as.numeric(mean_dat)), ylim = range(as.numeric(res_pred)), xlab = "True mean value", ylab = "Predicted value",
     asp = T, main = "Using our estimator")
lines(c(-100,100), rep(0,2), col = "red", lwd = 2)
lines(rep(0,2), c(-100,100), col = "red", lwd = 2)
lines(c(-100,100), c(-100,100), col = "red", lwd = 2)
points(as.numeric(mean_dat), as.numeric(res_pred), pch = 16, col = rgb(0,0,0,0.1))
lines(c(-100,100), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-100,100), col = "red", lty = 2)
lines(c(-100,100), c(-100,100), col = "red", lty = 2)

plot(NA, xlim = range(as.numeric(mean_dat)), ylim = range(as.numeric(res_pred_naive)), xlab = "True mean value", ylab = "Predicted value",
     asp = T, main = "Using naive estimator")
lines(c(-100,100), rep(0,2), col = "red", lwd = 2)
lines(rep(0,2), c(-100,100), col = "red", lwd = 2)
lines(c(-100,100), c(-100,100), col = "red", lwd = 2)
points(as.numeric(mean_dat), as.numeric(res_pred_naive), pch = 16, col = rgb(0,0,0,0.1))
lines(c(-100,100), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-100,100), col = "red", lty = 2)
lines(c(-100,100), c(-100,100), col = "red", lty = 2)

graphics.off()

####################

