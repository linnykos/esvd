rm(list=ls())
load("../../SOUP/data/zeisel.rda")

dat <- zeisel$counts

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 1000*log(dat + 1)

sparsity_vec <- apply(dat, 2, function(x){
  length(which(x!=0))
})

sum_vec <- colSums(dat)

quantile_sum <- quantile(sum_vec, prob = 0.95)

idx <- which(sum_vec <= quantile_sum)
dat <- dat[,idx]

idx <- which(colnames(dat) %in% zeisel$select.genes)
dat <- dat[,idx]

####################

# plot the matrix after we cluster
res_svd <- svd(dat)
u_mat <- res_svd$u[,1:8] %*% diag(res_svd$d[1:8])
v_mat <- res_svd$v[,1:8] %*% diag(res_svd$d[1:8])

u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 6, iter.max = 100, nstart = 10)
v_clust <- kmeans(v_mat_spherical, centers = 6, iter.max = 100, nstart = 10)

#####################

#reshuffle dat
row_idx <- unlist(lapply(1:6, function(x){
  sort(which(u_clust$cluster == x))
}))
col_idx <- unlist(lapply(1:6, function(x){
  sort(which(v_clust$cluster == x))
}))

dat <- dat[row_idx, col_idx]
dat2 <- dat

#####

# further rearrange based on quantiles
row_idx <- sapply(1:6, function(x){
  length(which(u_clust$cluster <= x))
})
col_idx <- sapply(1:6, function(x){
  length(which(v_clust$cluster <= x))
})
row_idx <- c(0, row_idx)
col_idx <- c(0, col_idx)

tmp <- as.numeric(dat)
tmp <- tmp[tmp!=0]
mid_val <- quantile(tmp, probs = 0.5)

# do columns first
for(i in 1:6){
  tmp <- dat[,c((col_idx[i]+1): col_idx[i+1])]

  low_percentage <- apply(tmp, 2, function(x){
    len <- length(x)
    x <- x[x!=0]
    length(which(x <= mid_val))/len
  })

  high_percentage <- apply(tmp, 2, function(x){
    len <- length(x)
    x <- x[x!=0]
    length(which(x >= mid_val))/len
  })

  high_low <- apply(tmp, 2, function(x){
    x <- x[x!=0]
    if(length(which(x >= mid_val)) > length(which(x < mid_val))) 1 else 0
  })

  num_val <- diff(range(col_idx[i]:col_idx[i+1]))
  stat_mat <- cbind(1:num_val, low_percentage, high_percentage, high_low)
  order_vec <- rep(0, num_val)
  bool <- TRUE
  counter_vec <- unique(unlist(lapply(1:ceiling(num_val/2), function(x){
    c(x, num_val-x+1)
  })))
  increment <- 1
  low_ignore <- F; high_ignore <- F
  while(TRUE){
    if(bool) {
      idx <- which.max(stat_mat[,2])
      if((stat_mat[idx,4] == 1 | low_ignore) & !high_ignore) {idx <- which.min(stat_mat[,3]); low_ignore = T}
    } else {
      idx <- which.max(stat_mat[,3])
      if((stat_mat[idx,4] == 0 | high_ignore) & !low_ignore) {idx <- which.min(stat_mat[,2]); high_ignore = T}
    }
    order_vec[counter_vec[increment]] <- stat_mat[idx,1]

    #remove from stat_mat
    stat_mat <- stat_mat[-idx,,drop = F]
    increment <- increment + 1
    bool <- !bool

    if(increment > num_val) break()
  }

  tmp <- tmp[,order_vec]
  dat[,c((col_idx[i]+1) : col_idx[i+1])] <- tmp
}


# then do rows
for(i in 1:6){
  tmp <- dat[c((row_idx[i]+1): row_idx[i+1]),]

  low_percentage <- apply(tmp, 1, function(x){
    len <- length(x)
    x <- x[x!=0]
    length(which(x <= mid_val))/len
  })

  high_percentage <- apply(tmp, 1, function(x){
    len <- length(x)
    x <- x[x!=0]
    length(which(x >= mid_val))/len
  })

  high_low <- apply(tmp, 1, function(x){
    x <- x[x!=0]
    if(length(which(x >= mid_val)) > length(which(x < mid_val))) 1 else 0
  })

  num_val <- diff(range(row_idx[i]:row_idx[i+1]))
  stat_mat <- cbind(1:num_val, low_percentage, high_percentage, high_low)
  order_vec <- rep(0, num_val)
  bool <- TRUE
  counter_vec <- unique(unlist(lapply(1:ceiling(num_val/2), function(x){
    c(x, num_val-x+1)
  })))
  increment <- 1
  low_ignore <- F; high_ignore <- F
  while(TRUE){
    if(bool) {
      idx <- which.max(stat_mat[,2])
      if((stat_mat[idx,4] == 1 | low_ignore) & !high_ignore) {idx <- which.min(stat_mat[,3]); low_ignore = T}
    } else {
      idx <- which.max(stat_mat[,3])
      if((stat_mat[idx,4] == 0 | high_ignore) & !low_ignore) {idx <- which.min(stat_mat[,2]); high_ignore = T}
    }
    order_vec[counter_vec[increment]] <- stat_mat[idx,1]

    #remove from stat_mat
    stat_mat <- stat_mat[-idx,,drop = F]
    increment <- increment + 1
    bool <- !bool

    if(increment > num_val) break()
  }

  tmp <- tmp[order_vec,]
  dat[c((row_idx[i]+1) : row_idx[i+1]),] <- tmp
}



#### set the colors
# col_vec <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)

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

png("../figure/experiment/4_image_reshuffled.png", height = 3600, width = 3600, res = 300, units = "px")
zlim <- range(dat)
break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)
image(.rotate(dat), breaks = break_vec, col = col_vec, asp = num_row/num_col,
      axes = F)

#put lines
row_idx <- sapply(1:5, function(x){
  1-length(which(u_clust$cluster <= x))/length(u_clust$cluster)
})
col_idx <- sapply(1:5, function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

for(i in row_idx){
  lines(c(0,1), rep(i, 2), lwd = 2, lty = 2)
}
for(i in col_idx){
  lines(rep(i, 2), c(0,1), lwd = 2, lty = 2)
}

graphics.off()


png("../figure/experiment/4_image.png", height = 3600, width = 3600, res = 300, units = "px")
zlim <- range(dat2)
break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)
image(.rotate(dat2), breaks = break_vec, col = col_vec, asp = num_row/num_col,
      axes = F)

#put lines
row_idx <- sapply(1:5, function(x){
  1-length(which(u_clust$cluster <= x))/length(u_clust$cluster)
})
col_idx <- sapply(1:5, function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

for(i in row_idx){
  lines(c(0,1), rep(i, 2), lwd = 2, lty = 2)
}
for(i in col_idx){
  lines(rep(i, 2), c(0,1), lwd = 2, lty = 2)
}

graphics.off()

table(zeisel$cell.info[,2], u_clust$cluster)

##########################

# inspect the covariance matrices

row_idx <- sapply(1:6, function(x){
  length(which(u_clust$cluster <= x))
})
row_idx <- c(0, row_idx)

col_idx <- sapply(1:5, function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

min_val <- Inf
max_val <- -Inf

for(i in 1:6){
  cov_mat <- stats::cov(dat[(row_idx[i]+1):row_idx[i+1],])
  min_val <- min(min_val, min(cov_mat))
  max_val <- max(max_val, max(cov_mat))
}

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)

for(i in 1:6){
  cov_mat <- stats::cov(dat[(row_idx[i]+1):row_idx[i+1],])

  break_vec <- quantile(as.numeric(cov_mat), probs = seq(0, 1, length.out = 20))
  break_vec[1] <- min_val; break_vec[20] <- max_val

  png(paste0("../figure/experiment/4_covariance_", i, "_reshuffled.png"), height = 3600, width = 3600, res = 300, units = "px")
  image(.rotate(cov_mat), breaks = break_vec, col = col_vec2, axes = F)

  for(j in col_idx){
    lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
  }
  for(j in col_idx){
    lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
  }

  graphics.off()
}

for(i in 1:6){
  cov_mat <- stats::cov(dat2[(row_idx[i]+1):row_idx[i+1],])

  break_vec <- quantile(as.numeric(cov_mat), probs = seq(0, 1, length.out = 20))
  break_vec[1] <- min_val; break_vec[20] <- max_val

  png(paste0("../figure/experiment/4_covariance_", i, ".png"), height = 3600, width = 3600, res = 300, units = "px")
  image(.rotate(cov_mat), breaks = break_vec, col = col_vec2, axes = F)

  for(j in col_idx){
    lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
  }
  for(j in col_idx){
    lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
  }

  graphics.off()
}

# plot entire covariance matrix
cov_mat <- stats::cov(dat2)
break_vec <- quantile(as.numeric(cov_mat), probs = seq(0, 1, length.out = 20))
break_vec[1] <- min_val; break_vec[20] <- max_val

png(paste0("../figure/experiment/4_covariance_all.png"), height = 3600, width = 3600, res = 300, units = "px")
image(.rotate(cov_mat), breaks = break_vec, col = col_vec2, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
}

graphics.off()

##############################

# new way to visualize data based on specific SBMs
pair_mat <- cbind(c(1,2,3,4,5), c(3,6,4,2,1), c(3, 4, 6, 4, 3))
colnames(pair_mat) <- c("Cell", "Gene", "Num.cluster")

row_idx <- sapply(1:6, function(x){
  length(which(u_clust$cluster <= x))
})
row_idx <- c(0, row_idx)

col_idx <- sapply(1:6, function(x){
  length(which(v_clust$cluster <= x))
})
col_idx <- c(0, col_idx)

png(paste0("../figure/experiment/4_block_eigenvalues.png"), height = 1800, width = 1200, res = 300, units = "px")
par(mfrow = c(3,2), mar = c(0,0,0,0))
for(i in pair_mat[,2]){
  idx <- pair_mat[which(pair_mat[,2] == i),1]
  cov_block <- stats::cov(dat2[(row_idx[idx]+1):row_idx[idx+1], (col_idx[i]+1):col_idx[i+1]])

  eig <- eigen(cov_block)

  plot(eig$values[1:50], pch = 16)
}
graphics.off()

png(paste0("../figure/experiment/4_block_eigenvectors.png"), height = 1800, width = 1200, res = 300, units = "px")
par(mfrow = c(3,2), mar = c(0,0,0,0))
for(i in 1:nrow(pair_mat)){
  cell_idx <- pair_mat[i,1]
  gene_idx <- pair_mat[i,2]
  cov_block <- stats::cov(dat2[(row_idx[cell_idx]+1):row_idx[cell_idx+1],
                               (col_idx[gene_idx]+1):col_idx[gene_idx+1]])

  eig <- eigen(cov_block)
  mat <- eig$vectors[,1:pair_mat[i,3]]
  mat <- t(apply(mat, 1, function(x){x/.l2norm(x)}))

  plot(mat[,1], mat[,2], pch = 16)
}
graphics.off()

for(i in 1:nrow(pair_mat)){
  cell_idx <- pair_mat[i,1]
  gene_idx <- pair_mat[i,2]

  png(paste0("../figure/experiment/4_block_covariance_", cell_idx, ".png"), height = 1200, width = 2400, res = 300, units = "px")
  par(mfrow = c(1,2), mar = c(1,1,1,1))

  cov_block <- stats::cov(dat2[(row_idx[cell_idx]+1):row_idx[cell_idx+1],
                               (col_idx[gene_idx]+1):col_idx[gene_idx+1]])

  threshold <- quantile(abs(as.numeric(cov_block)), probs = 0.99)
  cov_block[which(cov_block > threshold)] <- threshold
  cov_block[which(cov_block < -threshold)] <- -threshold

  break_vec <- quantile(as.numeric(cov_block), probs = seq(0, 1, length.out = 20))
  image(.rotate(cov_block), breaks = break_vec, col = col_vec2, asp = T, axes = F)

  eig <- eigen(cov_block)
  mat <- eig$vectors[,1:pair_mat[i,3]]
  mat <- t(apply(mat, 1, function(x){if(all(abs(x) < 1e-6)) x else x/.l2norm(x)}))

  set.seed(10)
  clust <- kmeans(mat, centers = pair_mat[i,3])
  idx <- unlist(lapply(1:pair_mat[i,3], function(x){
    which(clust$cluster == x)
  }))

  tmp_idx <- sapply(1:c(pair_mat[i,3]-1), function(x){
    length(which(clust$cluster <= x))/length(clust$cluster)
  })

  cov_block <- cov_block[idx, idx]
  image(.rotate(cov_block), breaks = break_vec, col = col_vec2, asp = T, axes = F)

  for(j in tmp_idx){
    lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
  }
  for(j in tmp_idx){
    lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
  }

  graphics.off()
}

## what about correlation?

for(i in 1:nrow(pair_mat)){
  cell_idx <- pair_mat[i,1]
  gene_idx <- pair_mat[i,2]

  png(paste0("../figure/experiment/4_block_correlation_", cell_idx, ".png"), height = 1200, width = 2400, res = 300, units = "px")
  par(mfrow = c(1,2), mar = c(1,1,1,1))

  cov_block <- stats::cor(dat2[(row_idx[cell_idx]+1):row_idx[cell_idx+1],
                               (col_idx[gene_idx]+1):col_idx[gene_idx+1]])

  threshold <- quantile(abs(as.numeric(cov_block)), probs = 0.99)
  cov_block[which(cov_block > threshold)] <- threshold
  cov_block[which(cov_block < -threshold)] <- -threshold

  break_vec <- quantile(as.numeric(cov_block), probs = seq(0, 1, length.out = 20))
  image(.rotate(cov_block), breaks = break_vec, col = col_vec2, asp = T, axes = F)

  eig <- eigen(cov_block)
  mat <- eig$vectors[,1:pair_mat[i,3]]
  mat <- t(apply(mat, 1, function(x){if(all(abs(x) < 1e-6)) x else x/.l2norm(x)}))

  set.seed(10)
  clust <- kmeans(mat, centers = pair_mat[i,3])
  idx <- unlist(lapply(1:pair_mat[i,3], function(x){
    which(clust$cluster == x)
  }))

  tmp_idx <- sapply(1:c(pair_mat[i,3]-1), function(x){
    length(which(clust$cluster <= x))/length(clust$cluster)
  })

  cov_block <- cov_block[idx, idx]
  image(.rotate(cov_block), breaks = break_vec, col = col_vec2, asp = T, axes = F)

  for(j in tmp_idx){
    lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
  }
  for(j in tmp_idx){
    lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
  }

  graphics.off()
}

##########################

# # try another strategy where we discard genes that don't seem to belong anywhere
# vec <- unlist(lapply(1:nrow(pair_mat), function(i){
#   cell_idx <- pair_mat[i,1]
#   gene_idx <- pair_mat[i,2]
#
#   tmp <- stats::cov(dat2[(row_idx[cell_idx]+1):row_idx[cell_idx+1],
#                   (col_idx[gene_idx]+1):col_idx[gene_idx+1]])
#   tmp[upper.tri(tmp, diag = T)]
# }))
#
# min_val <- quantile(vec, probs = 0.005)
# max_val <- quantile(vec, probs = 0.995)


##########################

# visualize correlation between two genes





