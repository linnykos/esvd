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
