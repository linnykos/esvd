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

# plot the matrix after we cluster
res_svd <- svd(dat)
u_mat <- res_svd$u[,1:8] %*% diag(res_svd$d[1:8])
v_mat <- res_svd$v[,1:8] %*% diag(res_svd$d[1:8])

u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 6, iter.max = 100, nstart = 10)
v_clust <- kmeans(v_mat_spherical, centers = 6, iter.max = 100, nstart = 10)

cell_type <- as.numeric(as.factor(zeisel$cell.info[,2]))
row_idx <- unlist(lapply(1:max(cell_type), function(x){
  sort(which(cell_type == x))
}))
col_idx <- unlist(lapply(1:6, function(x){
  sort(which(v_clust$cluster == x))
}))

dat <- dat[row_idx, col_idx]

tmp <- dat[dat != 0]
break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)

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

png("../figure/experiment/4_zeisel_true.png", height = 3600, width = 3600, res = 300, units = "px")
image(.rotate(dat), breaks = break_vec, col = col_vec, asp = nrow(dat)/ncol(dat),
      axes = F)

#put lines
row_idx <- sapply(1:(max(cell_type)-1), function(x){
  1-length(which(cell_type <= x))/length(u_clust$cluster)
})
col_idx <- sapply(1:(max(v_clust$cluster)-1), function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

for(i in row_idx){
  lines(c(0,1), rep(i, 2), lwd = 2, lty = 2)
}
for(i in col_idx){
  lines(rep(i, 2), c(0,1), lwd = 2, lty = 2)
}

graphics.off()

########################

# reshuffle the rows and columns
dat2 <- dat
row_idx <- sapply(1:max(cell_type), function(x){
  length(which(cell_type <= x))
})
col_idx <- sapply(1:max(v_clust$cluster), function(x){
  length(which(v_clust$cluster <= x))
})
row_idx <- c(0, row_idx)
col_idx <- c(0, col_idx)

tmp <- as.numeric(dat2)
tmp <- tmp[tmp!=0]
mid_val <- quantile(tmp, probs = 0.5)

# do columns first
for(i in 1:max(v_clust$cluster)){
  tmp <- dat2[,c((col_idx[i]+1): col_idx[i+1])]

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
  dat2[,c((col_idx[i]+1) : col_idx[i+1])] <- tmp
}


# then do rows
for(i in 1:max(cell_type)){
  tmp <- dat2[c((row_idx[i]+1): row_idx[i+1]),]

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
  dat2[c((row_idx[i]+1) : row_idx[i+1]),] <- tmp
}

png("../figure/experiment/4_zeisel_true_reshuffled.png", height = 3600, width = 3600, res = 300, units = "px")
image(.rotate(dat2), breaks = break_vec, col = col_vec, asp = nrow(dat)/ncol(dat),
      axes = F)

#put lines
row_idx <- sapply(1:(max(cell_type)-1), function(x){
  1-length(which(cell_type <= x))/length(u_clust$cluster)
})
col_idx <- sapply(1:(max(v_clust$cluster)-1), function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

for(i in row_idx){
  lines(c(0,1), rep(i, 2), lwd = 2, lty = 2)
}
for(i in col_idx){
  lines(rep(i, 2), c(0,1), lwd = 2, lty = 2)
}

graphics.off()


