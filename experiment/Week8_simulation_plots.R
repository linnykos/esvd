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

length(which(setting == 1))/length(which(setting >= 1))

##########

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

png("../figure/experiment/8_simulated_data.png", height = 1300, width = 2400, res = 300, units = "px")
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
