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

#########################

# compute the non-zero correlation matrix
d <- ncol(dat2)
cor_mat <- matrix(1, d, d)

dat_tmp <- dat2

nonzero_covariance <- function(i){
  if(i %% floor(ncol(combn_mat)/10) == 0) print('*')

  mat <- dat_tmp[,combn_mat[,i]]

  bool <- apply(mat, 1, function(x){(x[1] == 0 & x[2] != 0) || (x[1] != 0 & x[2] == 0)})
  if(all(bool)) return(0)

  mat <- mat[-which(bool),, drop = F]
  cor(mat[,1], mat[,2])
}

combn_mat <- combn(d, 2)
combn_mat <- cbind(combn_mat, rbind(1:d, 1:d))

doMC::registerDoMC(cores = 10)

cor_vec <- unlist(foreach::"%dopar"(foreach::foreach(i = 1:nrow(combn_mat), nonzero_covariance(i))))

save.image("../experiment/Week4_nonzero_covariance.RData")

########################

# do cell-type specific ones

cor_vec_list <- vector("list", 6)

row_idx <- sapply(1:6, function(x){
  length(which(u_clust$cluster <= x))
})
row_idx <- c(0, row_idx)

for(i in 1:6){
  dat_tmp <- dat2[(row_idx[idx]+1):row_idx[idx+1],]

  cor_vec_list[[i]] <- unlist(foreach::"%dopar"(foreach::foreach(i = 1:nrow(combn_mat), nonzero_covariance(i))))
}

save.image("../experiment/Week4_nonzero_covariance.RData")

