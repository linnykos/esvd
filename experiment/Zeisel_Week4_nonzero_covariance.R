rm(list=ls())

library(foreach)
library(doMC)

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

res_svd <- svd(dat)
u_mat <- res_svd$u[,1:8] %*% diag(res_svd$d[1:8])
v_mat <- res_svd$v[,1:8] %*% diag(res_svd$d[1:8])

.l2norm <- function(x){sqrt(sum(x^2))}

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

#########################

# compute the non-zero correlation matrix
d <- ncol(dat2)
n <- nrow(dat2)

doMC::registerDoMC(cores = 14)

nonzero_covariance <- function(i, row = F){
  print(i)
  if(row) mat <- t(dat2[combn_mat[,i],]) else mat <- dat2[,combn_mat[,i]]

  bool <- apply(mat, 1, function(x){(x[1] == 0 & x[2] != 0) || (x[1] != 0 & x[2] == 0)})
  if(all(bool)) return(0)

  if(any(bool)) mat <- mat[-which(bool),, drop = F]
  if(any(colSums(abs(mat)) == 0)) return(0)

  cov(mat[,1], mat[,2])
}

combn_mat <- combn(n, 2)
combn_mat <- cbind(combn_mat, rbind(1:n, 1:n))

i <- 1
cov_vec_cell <- unlist(foreach::"%dopar%"(foreach::foreach(i = 1:ncol(combn_mat)), nonzero_covariance(i, T)))

save(cov_vec_cell, file = "../experiment/Week4_nonzero_covariance_cell.RData")

################

combn_mat <- combn(d, 2)
combn_mat <- cbind(combn_mat, rbind(1:d, 1:d))

cov_vec_gene <- unlist(foreach::"%dopar%"(foreach::foreach(i = 1:ncol(combn_mat)), nonzero_covariance(i, F)))

save(cov_vec_gene, file = "../experiment/Week4_nonzero_covariance_gene.RData")
