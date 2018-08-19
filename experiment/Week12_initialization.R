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

gram_dat <- u_dat %*% t(v_dat)
n <- nrow(mean_dat)
d <- ncol(mean_dat)
dat <- matrix(0, n, d)
setting <- matrix(0, n, d)
nodropout_dat <- matrix(0, n, d)

prob_vec <- sample(seq(0.1, 0.45, length.out = ncol(gram_dat)))

set.seed(10)
for(i in 1:n){
  for(j in 1:d){
    val <- rnbinom(1, size = exp(gram_dat[i,j]), prob = prob_vec[j])
    bool <- 1
    setting[i,j] <- bool + 1
    dat[i,j] <- bool*val
    nodropout_dat[i,j] <- val
  }
}

mean_dat <- matrix(0, n, d)
var_dat <- matrix(0, n, d)

for(i in 1:n){
  for(j in 1:d){
    mean_dat[i,j] <- (1-prob_vec[j])*exp(gram_dat[i,j])/prob_vec[j]
    var_dat[i,j] <- (1-prob_vec[j])*exp(gram_dat[i,j])/prob_vec[j]^2
  }
}

k <- 3
tmp <- u_dat %*% t(v_dat)
res_svd <- svd(tmp)
u_dat_rescaled <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_dat_rescaled <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

length(which(dat == 0))/prod(dim(dat))

#####

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

png("../figure/experiment/12_simulated_data.png", height = 1300, width = 2400, res = 300, units = "px")
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

##########

# try initialization
# first, svd to predict the mean
res_svd <- svd(dat)
k <- 3
pred_mean <- res_svd$u[,1:k] %*% diag(res_svd$d[1:k]) %*% t(res_svd$v[,1:k])

#from predicted mean, infer variance
pred_var <- (pred_mean - dat)^2

#from variance and mean, predict theta
pred_theta <- sapply(1:ncol(dat), function(x){
  vec1 <- pred_mean[,x]; vec2 <- pred_var[,x]
  tmp <- as.data.frame(cbind(vec1, vec2))
  res <- lm(vec1 ~ vec2 - 1, data = tmp)
  coef(res)
})

#from predicted theta, recover exp gram matrix
pred_egram <- matrix(0, n, d)
for(i in 1:n){
  for(j in 1:d){
    pred_egram[i,j] <- pred_theta[j]*pred_mean[i,j]/(1-pred_theta[j])
  }
}

#from exp gram matrix, do an SVD
#janky non-negative transformation first
pred_egram2 <- log(exp(pred_egram)+1)
res_svd <- svd(log(pred_egram2))
pred_u <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
pred_v <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

#######################

png("../figure/experiment/12_mean_variance.png", height = 1000, width = 1800, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,4,4,1))
plot(as.numeric(mean_dat), as.numeric(pred_mean), asp = T, pch = 16,
     col = rgb(0,0,0,0.2), xlab = "True mean", ylab = "Predicted mean",
     main = "Mean")
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)

plot(as.numeric(var_dat), as.numeric(pred_var), asp = T, pch = 16,
     col = rgb(0,0,0,0.2), xlab = "True variance", ylab = "Predicted variance",
     main = "Variance")
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)
graphics.off()

png("../figure/experiment/12_theta.png", height = 1000, width = 1000, res = 300, units = "px")
plot(prob_vec, pred_theta, asp = T, pch = 16,
     col = rgb(0,0,0,0.5), xlab = "True theta", ylab = "Predicted theta",
     main = "Probability of failure")
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)
graphics.off()

png("../figure/experiment/12_gram.png", height = 1000, width = 1000, res = 300, units = "px")
plot(as.numeric(exp(gram_dat)), as.numeric(pred_egram2), asp = T, pch = 16,
     col = rgb(0,0,0,0.5), xlab = "True exp-gram", ylab = "Predicted exp-gram",
     main = "Exponentiated gram matrix")
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = "red", lwd = 2, lty = 2)
graphics.off()

png("../figure/experiment/12_latent_true.png", height = 1000, width = 1800, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,4,4,1))
plot(u_dat_rescaled[,1], u_dat_rescaled[,2], asp = T, pch = 16,
     xlab = "U[,1]", ylab = "U[,2]", main = "True Cells",
     col = c(1:3)[u_label])
lines(c(-1e5,1e5), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-1e5,1e5), col = "red", lwd = 2, lty = 2)

plot(v_dat_rescaled[,1], v_dat_rescaled[,2], asp = T, pch = 16,
     xlab = "V[,1]", ylab = "V[,2]", main = "True Genes",
     col = c(1:3)[v_label])
lines(c(-1e5,1e5), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-1e5,1e5), col = "red", lwd = 2, lty = 2)
graphics.off()

png("../figure/experiment/12_latent_predicted.png", height = 1000, width = 1800, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4,4,4,1))
plot(pred_u[,1], pred_u[,2], asp = T, pch = 16,
     xlab = "U[,1]", ylab = "U[,2]", main = "Predicted Cells",
     col = c(1:3)[u_label])
lines(c(-1e5,1e5), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-1e5,1e5), col = "red", lwd = 2, lty = 2)

plot(pred_v[,1], pred_v[,2], asp = T, pch = 16,
     xlab = "V[,1]", ylab = "V[,2]", main = "Predicted Genes",
     col = c(1:3)[v_label])
lines(c(-1e5,1e5), rep(0,2), col = "red", lwd = 2, lty = 2)
lines(rep(0,2), c(-1e5,1e5), col = "red", lwd = 2, lty = 2)
graphics.off()
