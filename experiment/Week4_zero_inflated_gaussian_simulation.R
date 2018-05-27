rm(list=ls())

# set up parameters
adj <- matrix(c(-.1, 0, -.1,
                0.15,-.25,-.25,
                -.75,-.25,0), 3, 3, byrow = T)
tmp <- svd(adj)
u_center <- t(tmp$u %*% diag(sqrt(tmp$d)))
v_center <- t(tmp$v %*% diag(sqrt(tmp$d)))

t(u_center) %*% v_center

u_num <- c(40, 20, 160)
u_label <- unlist(lapply(1:3, function(x){rep(x, u_num[x])}))
v_num <- c(60, 80, 200)
v_label <- unlist(lapply(1:3, function(x){rep(x, v_num[x])}))

# generate matrices
set.seed(10)
u_sig <- 0.05
u_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = u_num[x], mu = u_center[,x], Sigma = u_sig*diag(3))
}))
v_sig <- 0.05
v_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = v_num[x], mu = v_center[,x], Sigma = v_sig*diag(3))
}))

mean_dat <- u_dat %*% t(v_dat)
n <- nrow(mean_dat)
d <- ncol(mean_dat)
dat <- matrix(0, n, d)

dropout = 1-0.3
set.seed(10)
for(i in 1:n){
  for(j in 1:d){
    if(mean_dat[i,j] <= 0) {
      dat[i,j] <- 0
    } else {
      dat[i,j] <- rexp(1, rate = 1/mean_dat[i,j])*rbinom(1, size = 1, prob = dropout)
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

png("../figure/experiment/4_simulated_data.png", height = 1300, width = 2400, res = 300, units = "px")
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

################

# plot both covariance matrices

cov_d <- cov(dat)

col_idx <- cumsum(v_num)[1:2]/d

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_d), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_simulated_gene_covariance_true.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_d), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
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
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}
graphics.off()

###################

# let's try clustering now

res_svd <- svd(dat)
plot(res_svd$d[1:50])

k <- 3
u_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

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
  lines(c(0,1), rep(i, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 6, lty = 2)
}
for(i in col_idx){
  lines(rep(i, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(i, 2), c(0,1), lwd = 6, lty = 2)
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
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
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
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}
graphics.off()

#######

# what if we had the true means?

res_svd <- svd(mean_dat)
k <- 3
u_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))
# u_mat_spherical <- u_mat
# v_mat_spherical <- v_mat

set.seed(10)
u_clust2 <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
v_clust2 <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)

table(u_label, u_clust2$cluster)
table(v_label, v_clust2$cluster)

# manually reorder the clusters
u_clust2$cluster[u_clust2$cluster == 3] <- 4
u_clust2$cluster[u_clust2$cluster == 2] <- 5
u_clust2$cluster[u_clust2$cluster == 1] <- 6
u_clust2$cluster <- u_clust2$cluster - 3
v_clust2$cluster[v_clust2$cluster == 3] <- 5
v_clust2$cluster[v_clust2$cluster == 2] <- 4
v_clust2$cluster[v_clust2$cluster == 1] <- 6
v_clust2$cluster <- v_clust2$cluster - 3

# plot the reclustered mean_dat
row_idx <- unlist(lapply(1:max(u_clust2$cluster), function(x){
  sort(which(u_clust2$cluster == x))
}))
col_idx <- unlist(lapply(1:max(u_clust2$cluster), function(x){
  sort(which(v_clust2$cluster == x))
}))

dat3 <- dat[row_idx, col_idx]

tmp <- as.numeric(dat2)
tmp <- tmp[tmp!=0]

break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)

png("../figure/experiment/4_simulated_data_datadriven_meandat.png", height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(dat3), breaks = break_vec, col = col_vec, asp = nrow(dat2)/ncol(dat2),
      axes = F)

#put lines
row_idx <- sapply(1:2, function(x){
  1-length(which(u_clust2$cluster <= x))/length(u_clust2$cluster)
})
col_idx <- sapply(1:2, function(x){
  length(which(v_clust2$cluster <= x))/length(v_clust2$cluster)
})

for(i in row_idx){
  lines(c(0,1), rep(i, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(i, 2), lwd = 6, lty = 2)
}
for(i in col_idx){
  lines(rep(i, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(i, 2), c(0,1), lwd = 6, lty = 2)
}

graphics.off()

############################

# reorder the covariance matrices

cov_d <- cov(dat3)

col_idx <- sapply(1:(max(v_clust2$cluster)-1), function(x){
  length(which(v_clust2$cluster <= x))/length(v_clust2$cluster)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_d), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_simulated_gene_covariance_datadriven_meandat.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_d), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}
graphics.off()

###

cov_n <- cov(t(dat3))

col_idx <- sapply(1:(max(u_clust2$cluster)-1), function(x){
  length(which(u_clust2$cluster <= x))/length(u_clust2$cluster)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_d), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_simulated_cell_covariance_datadriven_meandat.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_n), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}
graphics.off()

#######

# plot the mean_dat

break_vec <- quantile(as.numeric(mean_dat), probs = seq(0, 1, length.out = 20))

png("../figure/experiment/4_simulated_data_mean.png", height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(mean_dat), breaks = break_vec, col = col_vec2, asp = nrow(dat)/ncol(dat),
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

###########

## "true" covariance matrices
# plot both covariance matrices

cov_d <- cov(mean_dat)

col_idx <- cumsum(v_num)[1:2]/d

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_d), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_simulated_gene_mean_covariance.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_d), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}
graphics.off()

###

cov_n <- cov(t(mean_dat))

col_idx <- cumsum(u_num)[1:2]/n

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_d), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_simulated_cell_mean_covariance.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_n), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}
graphics.off()

###############

# compare the gene covariance matrices, true, mean_dat, and dat, using real clustering
u_idx <- lapply(1:max(u_label), function(x){
  which(u_label == x)
})
v_idx <- lapply(1:max(v_label), function(x){
  which(u_label == x)
})

true_cov <- t(v_center)%*%v_center

tmp_cov <- v_dat %*% t(v_dat)
mean_cov <- matrix(0, 3, 3)
for(i in 1:3){
  for(j in i:3){
    mean_cov[i,j] <- mean(tmp_cov[v_idx[[i]], v_idx[[j]]])
    mean_cov[j,i] <- mean_cov[i,j]
  }
}

tmp_cov <- v_mat %*% t(v_mat)
dat_cov <- matrix(0, 3, 3)
for(i in 1:3){
  for(j in i:3){
    dat_cov[i,j] <- mean(tmp_cov[v_idx[[i]], v_idx[[j]]])
    dat_cov[j,i] <- dat_cov[i,j]
  }
}

true_cov
mean_cov
dat_cov
##
# do the same with cells
true_cov <- t(u_center)%*%u_center

tmp_cov <- u_dat %*% t(u_dat)
mean_cov <- matrix(0, 3, 3)
for(i in 1:3){
  for(j in i:3){
    mean_cov[i,j] <- mean(tmp_cov[v_idx[[i]], v_idx[[j]]])
    mean_cov[j,i] <- mean_cov[i,j]
  }
}

tmp_cov <- u_mat %*% t(u_mat)
dat_cov <- matrix(0, 3, 3)
for(i in 1:3){
  for(j in i:3){
    dat_cov[i,j] <- mean(tmp_cov[v_idx[[i]], v_idx[[j]]])
    dat_cov[j,i] <- dat_cov[i,j]
  }
}

true_cov
mean_cov
dat_cov

#############

# plot the sparsity vs sum

sparsity_vec <- apply(dat, 2, function(x){
  length(which(x!=0))
})

sum_vec <- colSums(dat)
n <- nrow(dat)

png("../figure/experiment/4_simulated_sparsity.png", height = 1800, width = 1800, res = 300, units = "px")
par(mar = c(5,5,1,1))
col_vec <- rep(rgb(0,0,0,0.5), n)
plot(sparsity_vec, sum_vec, col = col_vec, pch = 16,
     xlab = "Sparsity (Number of non-zeroes)", ylab = "Summation of counts",
     cex.lab = 2)
graphics.off()

################

# plot the confusion matrices
col_ramp <- colorRampPalette(c("white", rgb(0.803, 0.156, 0.211)))(10)
col_ramp <- sapply(col_ramp, function(x){paste0(x, "BF")})

png("../figure/experiment/4_simulated_table.png", height = 900, width = 1800, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(1,4,4,0))
u_tab1 <- table(u_label, u_clust$cluster)
u_tab2 <- t(apply(u_tab1, 1, function(x){x/sum(x)}))
num_row <- length(unique(u_label))
num_col <- length(unique(u_clust$cluster))


y_tic <- seq(-.5, .5+num_row-1, by = 1)*1/(num_row-1)
x_tic <- seq(-.5, .5+num_col-1, by = 1)*1/(num_col-1)
image(.rotate(u_tab2), col = col_ramp, asp = num_row/num_col, axes = F,
      xlim = c(min(x_tic)-.05, max(x_tic)+.05),
      ylim = c(min(y_tic)-.05, max(y_tic)+.05), main = "Cell confusion matrix\n(Based on X)")
for(y in y_tic){
  lines(range(x_tic), rep(y, 2))
}
for(x in x_tic){
  lines(rep(x, 2), range(y_tic))
}
text(par("usr")[3] - 0.1, seq(0,1,length.out = num_row), adj = 1, labels = rev(paste0("True ", 1:num_row)), xpd = TRUE,
     cex = 0.9)
text(seq(0, 1, length.out = num_col), par("usr")[1], labels = paste0("Est. ", 1:num_col), xpd = TRUE,
     cex = 0.9)

x_vec <- seq(0,1,length.out = num_col)
y_vec <- seq(0,1,length.out = num_row)
for(i in 1:length(x_vec)){
  for(j in 1:length(y_vec)){
    text(x_vec[i], 1-y_vec[j], label = u_tab1[j,i])
  }
}

###

v_tab1 <- table(v_label, v_clust$cluster)
v_tab2 <- t(apply(v_tab1, 1, function(x){x/sum(x)}))
num_row <- length(unique(v_label))
num_col <- length(unique(v_clust$cluster))

y_tic <- seq(-.5, .5+num_row-1, by = 1)*1/(num_row-1)
x_tic <- seq(-.5, .5+num_col-1, by = 1)*1/(num_col-1)
image(.rotate(v_tab2), col = col_ramp, asp = num_row/num_col, axes = F,
      xlim = c(min(x_tic)-.05, max(x_tic)+.05),
      ylim = c(min(y_tic)-.05, max(y_tic)+.05), main = "Gene confusion matrix\n(Based on X)")
for(y in y_tic){
  lines(range(x_tic), rep(y, 2))
}
for(x in x_tic){
  lines(rep(x, 2), range(y_tic))
}
text(par("usr")[3] - 0.1, seq(0,1,length.out = num_row), adj = 1, labels = rev(paste0("True ", 1:num_row)), xpd = TRUE,
     cex = 0.9)
text(seq(0, 1, length.out = num_col), par("usr")[1], labels = paste0("Est. ", 1:num_col), xpd = TRUE,
     cex = 0.9)

x_vec <- seq(0,1,length.out = num_col)
y_vec <- seq(0,1,length.out = num_row)
for(i in 1:length(x_vec)){
  for(j in 1:length(y_vec)){
    text(x_vec[i], 1-y_vec[j], label = v_tab1[j,i])
  }
}


graphics.off()

## now do the same for mean_dat


png("../figure/experiment/4_simulated_table_meandat.png", height = 900, width = 1800, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(1,4,4,0))
u_tab1 <- table(u_label, u_clust2$cluster)
u_tab2 <- t(apply(u_tab1, 1, function(x){x/sum(x)}))
num_row <- length(unique(u_label))
num_col <- length(unique(u_clust2$cluster))


y_tic <- seq(-.5, .5+num_row-1, by = 1)*1/(num_row-1)
x_tic <- seq(-.5, .5+num_col-1, by = 1)*1/(num_col-1)
image(.rotate(u_tab2), col = col_ramp, asp = num_row/num_col, axes = F,
      xlim = c(min(x_tic)-.05, max(x_tic)+.05),
      ylim = c(min(y_tic)-.05, max(y_tic)+.05), main = "Cell confusion matrix\n(Based on M)")
for(y in y_tic){
  lines(range(x_tic), rep(y, 2))
}
for(x in x_tic){
  lines(rep(x, 2), range(y_tic))
}
text(par("usr")[3] - 0.1, seq(0,1,length.out = num_row), adj = 1, labels = rev(paste0("True ", 1:num_row)), xpd = TRUE,
     cex = 0.9)
text(seq(0, 1, length.out = num_col), par("usr")[1], labels = paste0("Est. ", 1:num_col), xpd = TRUE,
     cex = 0.9)

x_vec <- seq(0,1,length.out = num_col)
y_vec <- seq(0,1,length.out = num_row)
for(i in 1:length(x_vec)){
  for(j in 1:length(y_vec)){
    text(x_vec[i], 1-y_vec[j], label = u_tab1[j,i])
  }
}

###

v_tab1 <- table(v_label, v_clust2$cluster)
v_tab2 <- t(apply(v_tab1, 1, function(x){x/sum(x)}))
num_row <- length(unique(v_label))
num_col <- length(unique(v_clust2$cluster))

y_tic <- seq(-.5, .5+num_row-1, by = 1)*1/(num_row-1)
x_tic <- seq(-.5, .5+num_col-1, by = 1)*1/(num_col-1)
image(.rotate(v_tab2), col = col_ramp, asp = num_row/num_col, axes = F,
      xlim = c(min(x_tic)-.05, max(x_tic)+.05),
      ylim = c(min(y_tic)-.05, max(y_tic)+.05), main = "Gene confusion matrix\n(Based on M)")
for(y in y_tic){
  lines(range(x_tic), rep(y, 2))
}
for(x in x_tic){
  lines(rep(x, 2), range(y_tic))
}
text(par("usr")[3] - 0.1, seq(0,1,length.out = num_row), adj = 1, labels = rev(paste0("True ", 1:num_row)), xpd = TRUE,
     cex = 0.9)
text(seq(0, 1, length.out = num_col), par("usr")[1], labels = paste0("Est. ", 1:num_col), xpd = TRUE,
     cex = 0.9)

x_vec <- seq(0,1,length.out = num_col)
y_vec <- seq(0,1,length.out = num_row)
for(i in 1:length(x_vec)){
  for(j in 1:length(y_vec)){
    text(x_vec[i], 1-y_vec[j], label = v_tab1[j,i])
  }
}


graphics.off()

