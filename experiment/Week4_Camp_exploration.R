load("../../SOUP/data/camp.rda")

dat <- camp$counts

dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 1000*log(dat + 1)

sparsity_vec <- apply(dat, 2, function(x){
  length(which(x!=0))
})

sum_vec <- colSums(dat)

quantile_sum <- quantile(sum_vec, prob = 0.95)

idx <- which(sum_vec <= quantile_sum)
dat <- dat[,idx]

idx <- which(colnames(dat) %in% camp$select.genes)
dat <- dat[,idx]

dat <- camp$counts[,which(colnames(camp$counts) %in% camp$select.genes)]
dim(dat)

##################

res_svd <- svd(dat)

plot(res_svd$d)
k <- 4
u_mat <- res_svd$u[,1:k] %*% diag(res_svd$d[1:k])
v_mat <- res_svd$v[,1:k] %*% diag(res_svd$d[1:k])

plot(u_mat[,1], u_mat[,2], pch = 16, asp = T)
plot(v_mat[,1], v_mat[,2], pch = 16, asp = T)

u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

# # try many number of clusters
# k_vec <- 2:20
# u_clust_vec <- sapply(k_vec, function(x){
#   set.seed(10)
#   kmeans(u_mat_spherical, centers = x, iter.max = 100, nstart = 10)$tot.withinss
# })
# v_clust_vec <- sapply(k_vec, function(x){
#   set.seed(10)
#   kmeans(v_mat_spherical, centers = x, iter.max = 100, nstart = 10)$tot.withinss
# })
#
# par(mfrow = c(1,2))
#
# plot(k_vec, u_clust_vec,  pch = 16, xlab = "K", ylab = "Total within ss", main = "Clustering for cells")
# plot(k_vec, v_clust_vec,pch = 16, xlab = "K", ylab = "Total within ss", main = "Clustering for genes")

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)

table(camp$cell.info[,2], u_clust$cluster)

#######################

# plot the data

#reshuffle dat
row_idx <- unlist(lapply(1:max(u_clust$cluster), function(x){
  sort(which(u_clust$cluster == x))
}))
col_idx <- unlist(lapply(1:max(u_clust$cluster), function(x){
  sort(which(v_clust$cluster == x))
}))

dat2 <- dat[row_idx, col_idx]

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

tmp <- as.numeric(dat2)
tmp <- tmp[tmp!=0]

break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)

png("../figure/experiment/4_camp_data.png", height = 2400, width = 2400, res = 300, units = "px")
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

##################

# do this with the real clusters

cell_type <- as.numeric(as.factor(camp$cell.info[,2]))
row_idx <- unlist(lapply(1:max(cell_type), function(x){
  which(cell_type == x)
}))
col_idx <- unlist(lapply(1:max(u_clust$cluster), function(x){
  sort(which(v_clust$cluster == x))
}))

dat3 <- dat[row_idx, col_idx]

tmp <- as.numeric(dat3)
tmp <- tmp[tmp!=0]

break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)

png("../figure/experiment/4_camp_data_sorted.png", height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(dat3), breaks = break_vec, col = col_vec, asp = nrow(dat3)/ncol(dat3),
      axes = F)

#put lines
row_idx <- sapply(1:(max(cell_type)-1), function(x){
  1-length(which(cell_type <= x))/length(cell_type)
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

###########


# let's visualize this covariance matrix
cov_mat <- stats::cov(t(dat2))

row_idx <- sapply(1:(max(u_clust$cluster)-1), function(x){
  length(which(u_clust$cluster <= x))/length(u_clust$cluster)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_mat), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_camp_cell_covariance.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_mat), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in row_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
}
for(j in row_idx){
  lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
}
graphics.off()

# and for genes

cov_mat <- stats::cov(dat2)

col_idx <- sapply(1:(max(v_clust$cluster)-1), function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_mat), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_camp_gene_covariance.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_mat), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
}
graphics.off()

#####################

tolerance <- quantile(dat[dat != 0], probs = 0.2)

nonzero_covariance <- function(i, tolerance = 0, row = F){
  print(i)
  if(row) mat <- t(dat[combn_mat[,i],]) else mat <- dat[,combn_mat[,i]]

  bool <- apply(mat, 1, function(x){(x[1] <= tolerance & x[2] > tolerance) ||
      (x[1] > tolerance & x[2] <= tolerance)})
  if(all(bool)) return(0)

  if(any(bool)) mat <- mat[-which(bool),, drop = F]
  if(any(colSums(abs(mat)) == 0)) return(0)

  cor(mat[,1], mat[,2])
}

d <- ncol(dat)
combn_mat <- combn(d, 2)
combn_mat <- cbind(combn_mat, rbind(1:d, 1:d))

doMC::registerDoMC(cores = 2)
i <- 1
cov_vec_gene <- unlist(foreach::"%dopar%"(foreach::foreach(i = 1:ncol(combn_mat)), nonzero_covariance(i, 20, F)))

cov_d <- matrix(0, d, d)
for(i in 1:ncol(combn_mat)){
  if(i %% floor(ncol(combn_mat)/10) == 0) print('*')

  cov_d[combn_mat[1,i], combn_mat[2,i]] <- cov_vec_gene[i]
  cov_d[combn_mat[2,i], combn_mat[1,i]] <- cov_vec_gene[i]
}

eig <- eigen(cov_d)
plot(eig$values[1:50])
k <- 2
v_mat <- eig$vectors[,1:k] %*% diag(eig$values[1:k])
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

k_vec <- 2:20
v_clust_vec <- sapply(k_vec, function(x){
  set.seed(10)
  kmeans(v_mat_spherical, centers = x, iter.max = 100, nstart = 10)$tot.withinss
})
plot(k_vec, v_clust_vec,pch = 16, xlab = "K", ylab = "Total within ss", main = "Clustering for genes")

set.seed(10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)

col_idx <- unlist(lapply(1:max(v_clust$cluster), function(x){
  sort(which(v_clust$cluster == x))
}))

cov_d <- cov_d[col_idx, col_idx]

col_idx <- sapply(1:(max(v_clust$cluster)-1), function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})
col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_d), probs = seq(0, 1, length.out = 20))
png(paste0("../figure/experiment/4_camp_gene_covariance_nonzero.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_d), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
}
graphics.off()

#############

# do the same but for cells

n <- nrow(dat)
combn_mat <- combn(n, 2)
combn_mat <- cbind(combn_mat, rbind(1:n, 1:n))

doMC::registerDoMC(cores = 2)
i <- 1
cov_vec_cell <- unlist(foreach::"%dopar%"(foreach::foreach(i = 1:ncol(combn_mat)), nonzero_covariance(i, 20, T)))

cov_n <- matrix(0, n, n)
for(i in 1:ncol(combn_mat)){
  if(i %% floor(ncol(combn_mat)/10) == 0) print('*')

  cov_n[combn_mat[1,i], combn_mat[2,i]] <- cov_vec_cell[i]
  cov_n[combn_mat[2,i], combn_mat[1,i]] <- cov_vec_cell[i]
}

eig <- eigen(cov_n)
plot(eig$values[1:50])

k <- 2
v_mat <- eig$vectors[,1:k] %*% diag(eig$values[1:k])
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

k_vec <- 2:20
v_clust_vec <- sapply(k_vec, function(x){
  set.seed(10)
  kmeans(v_mat_spherical, centers = x, iter.max = 100, nstart = 10)$tot.withinss
})
plot(k_vec, v_clust_vec,pch = 16, xlab = "K", ylab = "Total within ss", main = "Clustering for genes")

set.seed(10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)

col_idx <- unlist(lapply(1:max(v_clust$cluster), function(x){
  sort(which(v_clust$cluster == x))
}))

cov_n <- cov_n[col_idx, col_idx]

col_idx <- sapply(1:(max(v_clust$cluster)-1), function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})
col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_n), probs = seq(0, 1, length.out = 20))
png(paste0("../figure/experiment/4_camp_cell_covariance_nonzero.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_n), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in col_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
}
for(j in col_idx){
  lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
}
graphics.off()

#######################

for(i in 1:3){
  for(j in 1:3){
    idx1 <- which(u_clust$cluster == i)
    idx2 <- which(v_clust$cluster == j)

    tmp <- dat[idx1, idx2]
    tmp <- tmp[tmp != 0]
    hist(tmp, breaks = 50, col = "gray")
  }
}

#######################

#try nnfm

