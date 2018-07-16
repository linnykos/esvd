rm(list=ls())
adj <- matrix(c(10, 5, 2,
                15, 2, 5,
                5, 2, 10), 3, 3, byrow = T)
adj <- log(adj)
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
u_sig <- 0.01
u_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = u_num[x], mu = u_center[,x], Sigma = u_sig*diag(3))
}))
v_sig <- 0.01
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
dropout_coef <- c(-4, 1)
dropout_function <- dropout_create(dropout_coef[1], dropout_coef[2])

prob <- 0.7
set.seed(10)
for(i in 1:n){
  for(j in 1:d){
    val <- rnbinom(1, size = exp(mean_dat[i,j]), prob = prob)
    bool <- 1
    setting[i,j] <- bool + 1
    dat[i,j] <- bool*val
    nodropout_dat[i,j] <- val
  }
}

##################

library(zinbwave)
set.seed(10)
m <- zinbFit(dat, X = matrix(0, nrow = ncol(dat), ncol = 1),
             V = matrix(0, nrow = nrow(dat), ncol = 1), K = 3,
             verbose = T)
u_est <- m@W; v_est <- t(m@alpha_mu)

k <- 3
tmp <- u_dat %*% t(v_dat)
res_svd <- svd(tmp)
u_dat_rescaled <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_dat_rescaled <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

res_svd <- svd(dat)
k <- 3
u_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

######

png("../figure/experiment/10_cell_comparison_nodrop.png", height = 900, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(u_dat_rescaled[,1], u_dat_rescaled[,2], asp = T,
     pch = 16, col = c(1:3)[u_label], xlab = "U[,1]", ylab = "U[,2]",
     main = "Truth")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
plot(u_est[,1], u_est[,2], asp = T,
     pch = 16, col = c(1:3)[u_label], xlab = "U[,1]", ylab = "U[,1]",
     main = "Estimated via ZINB-WaVE")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
plot(u_mat[,1], u_mat[,2], asp = T,
     pch = 16, col = c(1:3)[u_label], xlab = "U[,1]", ylab = "U[,1]",
     main = "Estimated via SVD")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
graphics.off()


png("../figure/experiment/10_gene_comparison_nodrop.png", height = 900, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(v_dat_rescaled[,1], v_dat_rescaled[,2], asp = T,
     pch = 16, col = c(1:3)[v_label], xlab = "V[,1]", ylab = "V[,2]",
     main = "Truth")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
plot(v_est[,1], v_est[,2], asp = T,
     pch = 16, col = c(1:3)[u_label], xlab = "V[,1]", ylab = "V[,2]",
     main = "Estimated via ZINB-WaVE")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
plot(v_mat[,1], v_mat[,2], asp = T,
     pch = 16, col = c(1:3)[v_label], xlab = "V[,1]", ylab = "V[,1]",
     main = "Estimated via SVD")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
graphics.off()

nb_mean_dat <- exp(mean_dat)*(1-prob)/prob
nb_est_dat <- exp(u_est %*% t(v_est))*(1-prob)/prob
nb_est_svd_dat <- u_mat %*% t(v_mat)

png("../figure/experiment/10_mean_comparison_nodrop.png", height = 900, width = 2000, res = 300, units = "px")
par(mar = c(5,4,2,2), mfrow = c(1,2))
lim <- range(as.numeric(nb_mean_dat)); lim[1] <- -1
plot(as.numeric(nb_mean_dat), as.numeric(nb_est_dat), pch = 16, col = rgb(0,0,0,0.1),
     asp = T, xlab = "True mean", ylab = "Predicted mean",
     xlim = lim, ylim = lim, main = "Using ZINB-WaVE")
lines(c(-100,100), c(-100,100), col = "red", lwd = 2)

r_vec <- seq(0, max(mean_dat)*2, length.out = 100)
mean_vec <- exp(r_vec)*(1-prob)/prob
sd_vec <- sqrt(mean_vec/prob)
sd_plus <- mean_vec + sd_vec
sd_minus <- mean_vec - sd_vec
lines(mean_vec, sd_plus, col = "red", lwd = 2, lty = 2)
lines(mean_vec, sd_minus, col = "red", lwd = 2, lty = 2)

plot(as.numeric(nb_mean_dat), as.numeric(nb_est_svd_dat), pch = 16, col = rgb(0,0,0,0.1),
     asp = T, xlab = "True mean", ylab = "Predicted mean",
     xlim = lim, ylim = lim, main = "Using SVD")
lines(c(-100,100), c(-100,100), col = "red", lwd = 2)

r_vec <- seq(0, max(mean_dat)*2, length.out = 100)
mean_vec <- exp(r_vec)*(1-prob)/prob
sd_vec <- sqrt(mean_vec/prob)
sd_plus <- mean_vec + sd_vec
sd_minus <- mean_vec - sd_vec
lines(mean_vec, sd_plus, col = "red", lwd = 2, lty = 2)
lines(mean_vec, sd_minus, col = "red", lwd = 2, lty = 2)
graphics.off()

#############



