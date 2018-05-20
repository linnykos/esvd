rm(list=ls())

# set up parameters
set.seed(10)
tmp <- matrix(rnorm(9), 3, 3)
tmp <- tmp + t(tmp)
tmp <- eigen(tmp)$vectors
v_center <- matrix(0, 3, 3)
v_center[,1] <- tmp %*% c(-.5, .5, 1)
v_center[,2] <- tmp %*% c(1, 1, -.5)
v_center[,3] <- tmp %*% c(0, -.5, .5)

u_center <- matrix(0, 3, 3)
u_center[,1] <- solve(t(v_center), c(-0.25, -.1, -0.5))
u_center[,2] <- solve(t(v_center), c(-0.5, -.4, -0.25))
u_center[,3] <- solve(t(v_center), c(-0.75, -.4, -.5))

t(u_center)%*%v_center

u_num <- c(50, 80, 80)
u_label <- unlist(lapply(1:3, function(x){rep(x, u_num[x])}))
v_num <- c(80, 100, 200)
v_label <- unlist(lapply(1:3, function(x){rep(x, v_num[x])}))

# generate matrices
set.seed(10)
u_sig <- 0.1
u_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = u_num[x], mu = u_center[,x], Sigma = u_sig*diag(3))
}))
v_sig <- 0.1
v_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = v_num[x], mu = v_center[,x], Sigma = v_sig*diag(3))
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

png(paste0("../figure/experiment/4_simulated_gene_covariance_true.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_d), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
}
graphics.off()

###

cov_n <- cov(t(dat))

col_idx <- cumsum(u_num)[1:2]/n

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_d), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_simulated_cell_covariance_true.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_n), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
}
graphics.off()

###################

# let's try clustering now

res_svd <- svd(dat)
plot(res_svd$d[1:50])

k <- 3
u_mat <- res_svd$u[,1:k] %*% diag(res_svd$d[1:k])
v_mat <- res_svd$v[,1:k] %*% diag(res_svd$d[1:k])

plot(u_mat[,1], u_mat[,2], pch = 16, asp = T)
plot(v_mat[,1], v_mat[,2], pch = 16, asp = T)

u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)

table(u_label, u_clust$cluster)
table(v_label, v_clust$cluster)

#reshuffle dat
row_idx <- unlist(lapply(1:max(u_clust$cluster), function(x){
  sort(which(u_clust$cluster == x))
}))
col_idx <- unlist(lapply(1:max(u_clust$cluster), function(x){
  sort(which(v_clust$cluster == x))
}))

dat2 <- dat[row_idx, col_idx]

tmp <- as.numeric(dat2)
tmp <- tmp[tmp!=0]

break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)

png("../figure/experiment/4_simulated_data_datadriven.png", height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(dat2), breaks = break_vec, col = col_vec, asp = nrow(dat2)/ncol(dat2),
      axes = F)

#put lines
row_idx <- sapply(1:2, function(x){
  1-length(which(u_clust$cluster <= x))/length(u_clust$cluster)
})
col_idx <- sapply(1:2, function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

for(i in row_idx){
  lines(c(0,1), rep(i, 2), lwd = 2, lty = 2)
}
for(i in col_idx){
  lines(rep(i, 2), c(0,1), lwd = 2, lty = 2)
}

graphics.off()

############################

# reorder the covariance matrices

cov_d <- cov(dat2)

col_idx <- sapply(1:(max(v_clust$cluster)-1), function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_d), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_simulated_gene_covariance_datadriven.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_d), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
}
graphics.off()

###

cov_n <- cov(t(dat2))

col_idx <- sapply(1:(max(u_clust$cluster)-1), function(x){
  length(which(u_clust$cluster <= x))/length(u_clust$cluster)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_d), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_simulated_cell_covariance_datadriven.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_n), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
}
graphics.off()

#######

# what if we had the true means?

res_svd <- svd(mean_dat)
k <- 3
u_mat <- res_svd$u[,1:k] %*% diag(res_svd$d[1:k])
v_mat <- res_svd$v[,1:k] %*% diag(res_svd$d[1:k])

u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))
# u_mat_spherical <- u_mat
# v_mat_spherical <- v_mat

u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)

table(u_label, u_clust$cluster)
table(v_label, v_clust$cluster)
