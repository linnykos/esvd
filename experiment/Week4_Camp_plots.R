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

res_svd <- svd(dat)

k <- 4
u_mat <- res_svd$u[,1:k] %*% diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$v[,1:k] %*% diag(sqrt(res_svd$d[1:k]))

u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 3, iter.max = 100, nstart = 10)
v_clust <- kmeans(v_mat_spherical, centers = 3, iter.max = 100, nstart = 10)

##########################

# plot the data (group by major groups)
cell_type <- camp$cell.info[,2]
cell_type_coarse <- as.numeric(as.factor(sapply(cell_type, function(x){substr(x, 1, 1)})))

table(u_clust$cluster, cell_type_coarse)

row_idx <- unlist(lapply(1:max(cell_type_coarse), function(x){
  which(cell_type_coarse == x)
}))
col_idx <- unlist(lapply(1:max(v_clust$cluster), function(x){
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

tmp <- as.numeric(dat2)
tmp <- tmp[tmp!=0]

break_vec <- quantile(tmp, probs = seq(0, 1, length.out = 20))
break_vec <- c(-5, break_vec)

png("../figure/experiment/4_camp_data_final.png", height = 1300, width = 2400, res = 300, units = "px")
par(mar = rep(0.5,4))
image(.rotate(dat2), breaks = break_vec, col = col_vec, asp = nrow(dat2)/ncol(dat2),
      axes = F)

#put lines
row_idx <- sapply(1:(max(cell_type_coarse)-1), function(x){
  1-length(which(cell_type_coarse <= x))/length(cell_type)
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

############

# plot sparsity vs sum

sparsity_vec <- apply(dat, 2, function(x){
  length(which(x!=0))
})

sum_vec <- colSums(dat)
n <- nrow(dat)

png("../figure/experiment/4_camp_sparsity_final.png", height = 1500, width = 1500, res = 300, units = "px")
col_vec <- rep(rgb(0,0,0,0.5), n)
plot(sparsity_vec, sum_vec, col = col_vec, pch = 16,
     xlab = "Sparsity (Number of non-zeroes)", ylab = "Summation of counts")
graphics.off()

############

# plot the covariance matrices

cov_mat <- stats::cov(t(dat2))

row_idx <- sapply(1:(max(cell_type_coarse)-1), function(x){
  length(which(cell_type_coarse <= x))/length(cell_type_coarse)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_mat), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_camp_cell_covariance_final.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_mat), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in row_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 6, col = "white")
  lines(c(0,1), rep(1-j, 2), lwd = 6, lty = 2)
}
for(j in row_idx){
  lines(rep(j, 2), c(0,1), lwd = 6, col = "white")
  lines(rep(j, 2), c(0,1), lwd = 6, lty = 2)
}
graphics.off()

# and for genes

cov_mat <- stats::cov(dat2)

col_idx <- sapply(1:(max(v_clust$cluster)-1), function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_mat), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_camp_gene_covariance_final.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_mat), breaks = break_vec, col = col_vec2, asp = T, axes = F)

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

# we need to plot some cross sections. Let's do 4, two for genes, two for cells

# select gene pairs first
nonzero_covariance <- function(i, tolerance = 0, row = F){
  # print(i)
  if(row) mat <- t(dat[combn_mat[,i],]) else mat <- dat[,combn_mat[,i]]

  bool <- apply(mat, 1, function(x){(x[1] <= tolerance & x[2] > tolerance) ||
      (x[1] > tolerance & x[2] <= tolerance)})
  if(all(bool)) return(0)

  if(any(bool)) mat <- mat[-which(bool),, drop = F]
  if(any(colSums(abs(mat)) == 0)) return(0)

  cor(mat[,1], mat[,2])
}

nonzero_count <- function(i, tolerance = 0, row = F){
  if(row) mat <- t(dat[combn_mat[,i],]) else mat <- dat[,combn_mat[,i]]

  bool <- apply(mat, 1, function(x){all(x > tolerance)})

  length(which(bool))
}

set.seed(10)
d <- ncol(dat)
combn_mat <- combn(d, 2)
combn_mat <- combn_mat[,sample(1:ncol(combn_mat), 1000)]

cov_vec_gene_count <- sapply(1:ncol(combn_mat), function(x){nonzero_count(x, 0, F)})
cov_vec_gene_nonzero <- sapply(1:ncol(combn_mat), function(x){nonzero_covariance(x, 0, F)})
cov_vec_gene_full <- sapply(1:ncol(combn_mat), function(x){
  cor(dat[,combn_mat[1,x]], dat[,combn_mat[2,x]])
})

gene_pair1 <- combn_mat[,which.max(cov_vec_gene_count)]
tmp <- which(cov_vec_gene_count > quantile(cov_vec_gene_count, probs = 0.5))
gene_pair2 <- combn_mat[,tmp[which.max(abs(cov_vec_gene_nonzero[tmp] - cov_vec_gene_full[tmp]))]]

# next do cells
set.seed(10)
n <- nrow(dat)
combn_mat <- combn(n, 2)
combn_mat <- combn_mat[,sample(1:ncol(combn_mat), 1000)]

cov_vec_cell_count <- sapply(1:ncol(combn_mat), function(x){nonzero_count(x, 0, T)})
cov_vec_cell_nonzero <- sapply(1:ncol(combn_mat), function(x){nonzero_covariance(x, 0, T)})
cov_vec_cell_full <- sapply(1:ncol(combn_mat), function(x){
  cor(dat[combn_mat[1,x],], dat[combn_mat[2,x],])
})

cell_pair1 <- combn_mat[,which.max(cov_vec_cell_count)]
tmp <- which(cov_vec_cell_count > quantile(cov_vec_cell_count, probs = 0.5))
cell_pair2 <- combn_mat[,tmp[which.max(abs(cov_vec_cell_nonzero[tmp] - cov_vec_cell_full[tmp]))]]

# make the plot now

png(paste0("../figure/experiment/4_camp_scatterplot_final.png"), height = 1000, width = 3500, res = 300, units = "px")
par(mfrow = c(1,4), mar = c(5,5,4,1))

plot(dat[cell_pair1[1],], dat[cell_pair1[2],], pch = 16, col = rgb(0,0,0,0.5),
     xlab = paste0("Cell ", cell_pair1[1]), ylab = paste0("Cell ", cell_pair1[2]),
     main = paste0("Correlation: ", round(cor(dat[cell_pair1[1],], dat[cell_pair1[2],]), 2)),
     cex.lab = 2, cex.main = 2)

plot(dat[cell_pair2[1],], dat[cell_pair2[2],], pch = 16, col = rgb(0,0,0,0.5),
     xlab = paste0("Cell ", cell_pair2[1]), ylab = paste0("Cell ", cell_pair2[2]),
     main = paste0("Correlation: ", round(cor(dat[cell_pair2[1],], dat[cell_pair2[2],]), 2)),
     cex.lab = 2, cex.main = 2)

plot(dat[,gene_pair1[1]], dat[,gene_pair1[2]], pch = 16, col = rgb(0,0,0,0.5),
     xlab = paste0("Gene ", gene_pair1[1]), ylab = paste0("Gene ", gene_pair1[2]),
     main = paste0("Correlation: ", round(cor(dat[,gene_pair1[1]], dat[,gene_pair1[2]]), 2)),
     cex.lab = 2, cex.main = 2)

plot(dat[,gene_pair2[1]], dat[,gene_pair2[2]], pch = 16, col = rgb(0,0,0,0.5),
     xlab = paste0("Gene ", gene_pair2[1]), ylab = paste0("Gene ", gene_pair2[2]),
     main = paste0("Correlation: ", round(cor(dat[,gene_pair2[1]], dat[,gene_pair2[2]]), 2)),
     cex.lab = 2, cex.main = 2)

graphics.off()
