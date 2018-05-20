rm(list=ls())

# set up parameters
tmp <- matrix(rnorm(9), 3, 3)
tmp <- tmp + t(tmp)
u_mat <- eigen(tmp)$vectors


v_mat <- matrix(0, 3, 3)
v_mat[,1] <- solve(t(u_mat), c(-.2, -1.25, -1.5))
v_mat[,2] <- solve(t(u_mat), c(-.1, -1.25, -1.5))
v_mat[,3] <- solve(t(u_mat), c(-1, -.2, -.7))

t(u_mat)%*%v_mat

u_num <- c(50, 80, 80)
v_num <- c(80, 100, 200)

# generate matrices
set.seed(10)
sig <- 0.3
u_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = u_num[x], mu = u_mat[,x], Sigma = sig*diag(3))
}))
v_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = v_num[x], mu = v_mat[,x], Sigma = sig*diag(3))
}))

mean_dat <- u_dat %*% t(v_dat)
n <- nrow(mean_dat)
d <- ncol(mean_dat)
dat <- matrix(0, n, d)

set.seed(10)
for(i in 1:n){
  for(j in 1:d){
    if(mean_dat[i,j] <= 0) {
      dat[i,j] <- 0
    } else {
      dat[i,j] <- rexp(1, rate = 1/mean_dat[i,j])
    }
  }
}

################

# plot data

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

png("../figure/experiment/4_simulated_data.png", height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(dat), breaks = break_vec, col = col_vec, asp = nrow(dat)/ncol(dat),
      axes = F)

#put lines
row_idx <- 1-cumsum(u_num)[1:2]/n
col_idx <- cumsum(v_num)[1:2]/d

for(i in row_idx){
  lines(c(0,1), rep(i, 2), lwd = 2, lty = 2)
}
for(i in col_idx){
  lines(rep(i, 2), c(0,1), lwd = 2, lty = 2)
}

graphics.off()

################

# plot both covariance matrices

cov_d <- cov(dat)

col_idx <- cumsum(v_num)[1:2]/d

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_d), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_simulated_gene_covariance.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_d), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
}
graphics.off()
