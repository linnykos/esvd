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

#####################

#reshuffle dat
row_idx <- unlist(lapply(1:6, function(x){
  sort(which(u_clust$cluster == x))
}))
col_idx <- unlist(lapply(1:6, function(x){
  sort(which(v_clust$cluster == x))
}))

dat <- dat[row_idx, col_idx]

#reshuffle the labels as well
cell_type <- zeisel$cell.info[row_idx,2]

load("../experiment/Week4_nonzero_covariance_cell.RData")

###########

# form covariance matrix

n <- nrow(dat)
stopifnot(length(cov_vec_cell) == n*(n-1)/2 + n)

combn_mat <- combn(n, 2)
combn_mat <- cbind(combn_mat, rbind(1:n, 1:n))

cov_n <- matrix(0, n, n)
for(i in 1:ncol(combn_mat)){
  if(i %% floor(ncol(combn_mat)/10) == 0) print('*')

  cov_n[combn_mat[1,i], combn_mat[2,i]] <- cov_vec_cell[i]
  cov_n[combn_mat[2,i], combn_mat[1,i]] <- cov_vec_cell[i]
}

eig <- eigen(cov_n)

plot(eig$vectors[,1], eig$vectors[,2], pch = 16, col = rgb(0,0,0,0.2), asp = T)
plot(eig$values[1:50])

k <- 5
u_mat <- eig$vectors[,1:k] %*% diag(eig$values[1:k])
u_mat_spherical <- t(apply(u_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
u_clust <- kmeans(u_mat_spherical, centers = 6, iter.max = 100, nstart = 10)

.rotate <- function(mat){t(mat)[,nrow(mat):1]}
u_tab <- table(cell_type, u_clust$cluster)
u_tab2 <- t(apply(u_tab, 1, function(x){x/sum(x)}))


png("../figure/experiment/4_table_nonzero_cell.png", height = 1200, width = 2400, res = 300, units = "px")
par(mar = c(1,4,4,0))
col_ramp <- colorRampPalette(c("white", rgb(0.803, 0.156, 0.211)))(100)

y_tic <- seq(-.5, 6.5, by = 1)*1/6
x_tic <- seq(-.5, 5.5, by = 1)*1/5
image(.rotate(u_tab2), col = col_ramp, asp = 7/5, axes = F,
      xlim = c(min(x_tic)-.05, max(x_tic)+.05),
      ylim = c(min(y_tic)-.05, max(y_tic)+.05), main = "Cluster table for Cells\n(Count)")
for(y in y_tic){
  lines(range(x_tic), rep(y, 2))
}
for(x in x_tic){
  lines(rep(x, 2), range(y_tic))
}
text(par("usr")[3] - 0.1, seq(0,1,length.out = 7), adj = 1, labels = rev(rownames(u_tab)), xpd = TRUE,
     cex = 0.9)
text(seq(0, 1, length.out = 6), par("usr")[1] + 0.3, labels = paste0("C", 1:6), xpd = TRUE,
     cex = 0.9)

x_vec <- seq(0,1,length.out = 6)
y_vec <- seq(0,1,length.out = 7)
for(i in 1:length(x_vec)){
  for(j in 1:length(y_vec)){
    text(x_vec[i], 1-y_vec[j], label = u_tab[j,i])
  }
}
graphics.off()

##################

# let's try the same but for correlation
cor_n <- cov_n
cor_n <- diag(1/sqrt(diag(cor_n))) %*% cor_n %*% diag(1/sqrt(diag(cor_n)))
# ... why are some elements larger than 1?

eig <- eigen(cor_n)

k <- 2
u_mat <- eig$vectors[,1:k] %*% diag(eig$values[1:k])

col_idx <- as.numeric(as.factor(zeisel$cell.info[,2]))
plot(u_mat[,1], u_mat[,2], col = col_idx, pch = 16)

###################

# let's visualize this covariance matrix
idx <- unlist(lapply(1:max(u_clust$cluster), function(x){which(u_clust$cluster == x)}))
cov_n <- cov_n[idx,idx]

row_idx <- sapply(1:(max(u_clust$cluster)-1), function(x){
  length(which(u_clust$cluster <= x))/length(u_clust$cluster)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_n), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_cell_nonzero_cvoariance.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_n), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in row_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
}
for(j in row_idx){
  lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
}
graphics.off()

########################

# now do gene

load("../experiment/Week4_nonzero_covariance_gene.RData")

d <- ncol(dat)
stopifnot(length(cov_vec_gene) == d*(d-1)/2 + d)

combn_mat <- combn(d, 2)
combn_mat <- cbind(combn_mat, rbind(1:d, 1:d))

cov_d <- matrix(0, d, d)
for(i in 1:ncol(combn_mat)){
  if(i %% floor(ncol(combn_mat)/10) == 0) print('*')

  cov_d[combn_mat[1,i], combn_mat[2,i]] <- cov_vec_gene[i]
  cov_d[combn_mat[2,i], combn_mat[1,i]] <- cov_vec_gene[i]
}

eig <- eigen(cov_d)

plot(eig$vectors[,1], eig$vectors[,2], pch = 16, col = rgb(0,0,0,0.2), asp = T)
plot(eig$values[1:50])

k <- 6
v_mat <- eig$vectors[,1:k] %*% diag(eig$values[1:k])
v_mat_spherical <- t(apply(v_mat, 1, function(x){x/.l2norm(x)}))

set.seed(10)
v_clust <- kmeans(v_mat_spherical, centers = 6, iter.max = 100, nstart = 10)

# plotting

idx <- unlist(lapply(1:max(v_clust$cluster), function(x){which(v_clust$cluster == x)}))
cov_d <- cov_d[idx,idx]

row_idx <- sapply(1:(max(v_clust$cluster)-1), function(x){
  length(which(v_clust$cluster <= x))/length(v_clust$cluster)
})

col_vec2 <- colorRampPalette(c(rgb(0.584, 0.858, 0.564), rgb(0.803, 0.156, 0.211)))(19)
break_vec <- quantile(as.numeric(cov_d), probs = seq(0, 1, length.out = 20))

png(paste0("../figure/experiment/4_gene_nonzero_cvoariance.png"), height = 2400, width = 2400, res = 300, units = "px")
image(.rotate(cov_d), breaks = break_vec, col = col_vec2, asp = T, axes = F)

for(j in row_idx){
  lines(c(0,1), rep(1-j, 2), lwd = 2, lty = 2)
}
for(j in row_idx){
  lines(rep(j, 2), c(0,1), lwd = 2, lty = 2)
}
graphics.off()

## try gene correlation instead?
cor_d <- cov_d
cor_d <- diag(1/sqrt(diag(cor_d))) %*% cor_d %*% diag(1/sqrt(diag(cor_d)))

eig <- eigen(cor_d)

