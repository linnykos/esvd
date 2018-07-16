rm(list=ls())
adj <- matrix(c(10, 5, 2,
                5, 2, 1,
                2, 1, 5), 3, 3, byrow = T)
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
dropout_coef <- c(0, 3/4)
dropout_function <- dropout_create(dropout_coef[1], dropout_coef[2])

prob <- 0.7
set.seed(10)
for(i in 1:n){
  for(j in 1:d){
    val <- rnbinom(1, size = exp(mean_dat[i,j]), prob = prob)
    bool <- rbinom(1, 1, prob = dropout_function(val))
    setting[i,j] <- bool + 1
    dat[i,j] <- bool*val
    nodropout_dat[i,j] <- val
  }
}

length(intersect(which(setting == 2), which(dat == 0)))/prod(dim(dat))
length(which(setting == 1))/prod(dim(dat))
table(as.numeric(setting))
quantile(as.numeric(dat))
quantile(as.numeric(nodropout_dat))
length(which(nodropout_dat == 0))/prod(dim(dat))
length(intersect(which(nodropout_dat == 0), which(setting == 1)))/prod(dim(dat))
dropout_function(0)

x <- seq(0, max(nodropout_dat), length.out = 100)
y <- sapply(x, dropout_function)

x2 <- sapply(1:100, function(x){quantile(nodropout_dat, prob = x/100)})
y2 <- sapply(x2, dropout_function)

par(mfrow = c(1,2))
plot(x,y)
plot(jitter(x2), y2)

##################

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

colorRamp_custom <- function(vec1, vec2, length){
  mat <- matrix(0, nrow = length, ncol = 3)
  for(i in 1:3){
    mat[,i] <- seq(vec1[i], vec2[i], length.out = length)
  }

  luminosity_vec <- apply(mat, 1, function(x){
    0.2126*x[1] + 0.7152*x[2] + 0.0722*x[3]
  })

  target_luminosity <- mean(c(luminosity_vec[1], luminosity_vec[length]))

  mat <- t(sapply(1:nrow(mat), function(x){
    factor <- min(c(target_luminosity/luminosity_vec[x], 1/mat[x,]))
    mat[x,] * factor
  }))

  apply(mat, 1, function(x){
    rgb(x[1], x[2], x[3])
  })
}

col_vec <- colorRamp_custom(c(0.584, 0.858, 0.564), c(0.803, 0.156, 0.211), 19)
col_vec <- c("white", col_vec)

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

tmp <- as.numeric(dat)
tmp <- tmp[tmp!=0]

break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)

png("../figure/experiment/10_simulated_data.png", height = 1300, width = 2400, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(dat), breaks = break_vec, col = col_vec, asp = nrow(dat)/ncol(dat),
      axes = F)

#put lines
row_idx <- 1-cumsum(u_num)[1:2]/n
col_idx <- cumsum(v_num)[1:2]/d

for(i in row_idx){
  lines(c(0,1), rep(i, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 6, lty = 2)
}
for(i in col_idx){
  lines(rep(i, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(i, 2), c(0,1), lwd = 6, lty = 2)
}

graphics.off()

png("../figure/experiment/10_simulated_data_hist.png", height = 1100, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(5,4,2,0))
hist(as.numeric(dat), breaks = 20, col = "gray",
     main = "Histogram\nof observed data", xlab = "Value")
hist(as.numeric(nodropout_dat), breaks = 20, col = "gray",
     main = "Histogram\nof uncensored data", xlab = "Value")
graphics.off()


##############

# now try to estimate

library(zinbwave)
set.seed(10)
m <- zinbFit(dat, X = matrix(0, nrow = ncol(dat), ncol = 1),
             V = matrix(0, nrow = nrow(dat), ncol = 1), K = 3)
# m@W
# m@alpha_mu
# m@alpha_pi

######

k <- 3
tmp <- u_dat %*% t(v_dat)
res_svd <- svd(tmp)
u_dat_rescaled <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_dat_rescaled <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

# rescale the estimates by PCA
# tmp <- m@W %*% m@alpha_mu
# res_svd <- svd(tmp)
# u_est <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
# v_est <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
u_est <- m@W; v_est <- t(m@alpha_mu)


# try good-ol spectral clustering
res_svd <- svd(dat)
k <- 3
u_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
.l2norm <- function(x){sqrt(sum(x^2))}
u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))
set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)

###########################

png("../figure/experiment/10_true_latent.png", height = 1100, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(u_dat_rescaled[,1], u_dat_rescaled[,2], asp = T,
     pch = 16, col = c(1:3)[u_label], xlab = "U[,1]", ylab = "U[,2]",
     main = "Truth (Cells)")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)

plot(v_dat_rescaled[,1], v_dat_rescaled[,2], asp = T,
     pch = 16, col = c(1:3)[v_label], xlab = "V[,1]", ylab = "V[,2]",
     main = "Truth (Genes)")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
graphics.off()

#############################

# compare the means
nb_mean_dat <- exp(mean_dat)*(1-prob)/prob
nb_est_dat <- exp(u_est %*% t(v_est))*(1-prob)/prob
nb_est_svd_dat <- u_mat %*% t(v_mat)

# do an SVD on log(X), and then transform it back?
res_svd <- svd(log(dat + 1))
u_dat_log <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_dat_log <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
pred_mat_svd <- exp(u_dat_log %*% t(v_dat_log))

res_svd <- svd(nb_mean_dat)
u_dat_reparameterized <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_dat_reparameterized <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

res_svd <- svd(nb_est_dat)
u_est_reparameterized <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_est_reparameterized <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))


# side by side plots
# cells first

png("../figure/experiment/10_cell_comparison.png", height = 900, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(u_dat_reparameterized[,1], -u_dat_reparameterized[,2], asp = T,
     pch = 16, col = c(1:3)[u_label], xlab = "U[,1]", ylab = "U[,2]",
     main = "Truth\n(Reparameterized)")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
plot(u_est_reparameterized[,1], -u_est_reparameterized[,2], asp = T,
     pch = 16, col = c(1:3)[u_label], xlab = "U[,1]", ylab = "U[,1]",
     main = "Estimated via ZINB-WaVE\n(Reparameterized)")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
plot(u_mat[,1], u_mat[,2], asp = T,
     pch = 16, col = c(1:3)[u_label], xlab = "U[,1]", ylab = "U[,1]",
     main = "Estimated via SVD")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
graphics.off()

##

png("../figure/experiment/10_gene_comparison.png", height = 900, width = 2400, res = 300, units = "px")
par(mfrow = c(1,3))
plot(v_dat_reparameterized[,1], -v_dat_reparameterized[,2], asp = T,
     pch = 16, col = c(1:3)[v_label], xlab = "V[,1]", ylab = "V[,2]",
     main = "Truth\n(Reparameterized)")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
plot(v_est_reparameterized[,1], -v_est_reparameterized[,2], asp = T,
     pch = 16, col = c(1:3)[v_label], xlab = "V[,1]", ylab = "V[,2]",
     main = "Estimated via ZINB-WaVE\n(Reparameterized)")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
plot(v_mat[,1], v_mat[,2], asp = T,
     pch = 16, col = c(1:3)[v_label], xlab = "V[,1]", ylab = "V[,1]",
     main = "Estimated via SVD")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
graphics.off()

################

#par(mfrow = c(1,2))
#plot(as.numeric(nodropout_dat), as.numeric(nb_mean_dat), pch = 16, col = rgb(0,0,0,0.3))
png("../figure/experiment/10_mean_comparison.png", height = 900, width = 2000, res = 300, units = "px")
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
lines(c(-100, 100), rep(0, 2), col = "red", lwd = 2, lty = 2)

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
lines(c(-100, 100), rep(0, 2), col = "red", lwd = 2, lty = 2)
graphics.off()

#################


#par(mfrow = c(1,2))
#plot(as.numeric(nodropout_dat), as.numeric(nb_mean_dat), pch = 16, col = rgb(0,0,0,0.3))
png("../figure/experiment/10_observed_comparison.png", height = 900, width = 2000, res = 300, units = "px")
par(mar = c(5,4,2,2), mfrow = c(1,2))
lim <- range(as.numeric(dat)); lim[1] <- -1
plot(as.numeric(dat), as.numeric(nb_est_dat), pch = 16, col = rgb(0,0,0,0.1),
     asp = T, xlab = "True mean", ylab = "Predicted mean",
     xlim = lim, ylim = lim, main = "Using ZINB-WaVE")
lines(c(-100,100), c(-100,100), col = "red", lwd = 2)
lines(c(-100, 100), rep(0, 2), col = "red", lwd = 2, lty = 2)

plot(as.numeric(dat), as.numeric(nb_est_svd_dat), pch = 16, col = rgb(0,0,0,0.1),
     asp = T, xlab = "True mean", ylab = "Predicted mean",
     xlim = lim, ylim = lim, main = "Using SVD")
lines(c(-100,100), c(-100,100), col = "red", lwd = 2)
lines(c(-100, 100), rep(0, 2), col = "red", lwd = 2, lty = 2)
graphics.off()

##############

plot(sapply(as.numeric(nodropout_dat), dropout_function),
     as.numeric(getPi(m)), pch = 16, asp = T)

#################

png("../figure/experiment/10_cell_comparison2.png", height = 1000, width = 1800, res = 300, units = "px")
par(mfrow = c(1,2))
plot(u_dat_reparameterized[,1], -u_dat_reparameterized[,2], asp = T,
     pch = 16, col = c(1:3)[u_label], xlab = "U[,1]", ylab = "U[,2]",
     main = "Truth (Reparameterized)")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)

plot(u_mat[,1], u_mat[,2], asp = T,
     pch = 16, col = c(1:3)[u_label], xlab = "U[,1]", ylab = "U[,1]",
     main = "Estimated via SVD")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
graphics.off()

png("../figure/experiment/10_gene_comparison2.png", height = 1000, width = 1800, res = 300, units = "px")
par(mfrow = c(1,2))
plot(v_dat_reparameterized[,1], -v_dat_reparameterized[,2], asp = T,
     pch = 16, col = c(1:3)[v_label], xlab = "V[,1]", ylab = "V[,2]",
     main = "Truth (Reparameterized)")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)

plot(v_mat[,1], v_mat[,2], asp = T,
     pch = 16, col = c(1:3)[v_label], xlab = "V[,1]", ylab = "V[,1]",
     main = "Estimated via SVD")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
graphics.off()

###########################

# do a janky heurisitc: predictions above a certain threshold are viewed
# as a true dropout if the recorded value is 0. Let's see how well this works

range_vec <- range(nb_est_svd_dat)
thres_vec <- seq(0, range_vec[2], length.out = 100)
tpr_vec <- rep(NA, length(thres_vec))
fpr_vec <- rep(NA, length(thres_vec))

idx_true <- which(setting == 1)
idx_false <- which(setting == 2)
idx_zero <- which(dat == 0)
for(i in 1:length(thres_vec)){
  idx <- which(nb_est_svd_dat >= thres_vec[i])
  idx <- intersect(idx, idx_zero)
  tpr_vec[i] <- length(intersect(idx, idx_true))/length(idx_true)
  fpr_vec[i] <- length(intersect(idx, idx_false))/length(idx_false)
}

tpr_vec_log <- rep(NA, length(thres_vec))
fpr_vec_log <- rep(NA, length(thres_vec))

for(i in 1:length(thres_vec)){
  idx <- which(pred_mat_svd >= thres_vec[i])
  idx <- intersect(idx, idx_zero)
  tpr_vec_log[i] <- length(intersect(idx, idx_true))/length(idx_true)
  fpr_vec_log[i] <- length(intersect(idx, idx_false))/length(idx_false)
}


png("../figure/experiment/10_roc_dropout.png", height = 1100, width = 1800, res = 300, units = "px")
par(mar = c(5,4,4,1), mfrow = c(1,2))
plot(fpr_vec, tpr_vec, xlim = c(0,1), ylim = c(0,1), asp = T,
     pch = 16, xlab = "False positive rate", ylab = "True positive rate",
     main = "ROC for SVD")
lines(c(0,1), c(0,1), col = "red", lwd = 2, lty = 2)

plot(fpr_vec_log, tpr_vec_log, xlim = c(0,1), ylim = c(0,1), asp = T,
     pch = 16, xlab = "False positive rate", ylab = "True positive rate",
     main = "ROC for SVD of log")
lines(c(0,1), c(0,1), col = "red", lwd = 2, lty = 2)
graphics.off()

vec <- as.numeric(nodropout_dat)
idx <- order(vec)
plot(vec[idx], col = c(1,2)[as.numeric(setting)[idx]], pch = 16)

plot(tpr_vec)

##############################


png("../figure/experiment/10_log_analysis.png", height = 800, width = 2000, res = 300, units = "px")
par(mfrow = c(1,3))
plot(as.numeric(nb_mean_dat), as.numeric(pred_mat_svd), pch = 16, col = rgb(0,0,0,0.1),
     asp = T, xlab = "True mean", ylab = "Predicted mean",
     xlim = lim, ylim = lim, main = "Using SVD of log")

lines(c(-100,100), c(-100,100), col = "red", lwd = 2)

r_vec <- seq(0, max(mean_dat)*2, length.out = 100)
mean_vec <- exp(r_vec)*(1-prob)/prob
sd_vec <- sqrt(mean_vec/prob)
sd_plus <- mean_vec + sd_vec
sd_minus <- mean_vec - sd_vec
lines(mean_vec, sd_plus, col = "red", lwd = 2, lty = 2)
lines(mean_vec, sd_minus, col = "red", lwd = 2, lty = 2)
lines(c(-100, 100), rep(0, 2), col = "red", lwd = 2, lty = 2)

plot(u_dat_log[,1], u_dat_log[,2], asp = T,
     pch = 16, col = c(1:3)[u_label], xlab = "U[,1]", ylab = "U[,1]",
     main = "Estimated via SVD of log")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)

plot(v_dat_log[,1], v_dat_log[,2], asp = T,
     pch = 16, col = c(1:3)[v_label], xlab = "V[,1]", ylab = "V[,1]",
     main = "Estimated via SVD of log")
lines(c(-5,5), rep(0,2), col = "red", lty = 2)
lines(rep(0,2), c(-5,5), col = "red", lty = 2)
graphics.off()
