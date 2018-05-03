rm(list=ls())
set.seed(10)

k <- 25
d <- 100
n <- 50

mat1 <- sapply(1:k, function(x){
  rexp(n, rate = 10*1/x)
})
mat2 <- sapply(1:k, function(x){
  rexp(d, rate = 10*1/x)
})

mean_mat <- mat1 %*% t(mat2)
plot(as.numeric(mean_mat))

obs_mat <- mean_mat
for(i in 1:n){
  for(j in 1:d){
    obs_mat[i,j] <- rexp(1, rate = 1/mean_mat[i,j])
  }
}

##########################
# hypothetical lines
plot(as.numeric(mean_mat), as.numeric(obs_mat), asp = T)
mean_val <- seq(0, 500, by = 1)
lines(mean_val, mean_val, col = "red", lwd = 2)
lines(mean_val, 2*mean_val, col = "red", lwd = 2, lty = 2)
lines(mean_val, rep(0, length(mean_val)), col = "red", lwd = 2, lty = 2)
###########################

# compute hypothetical svd
mean_svd <- svd(mean_mat)
plot(mean_svd$d)

#########################

obs_svd <- svd(obs_mat)
plot(obs_svd$d)
pred_mean_mat <- obs_svd$u[,1,drop=F] %*% obs_svd$d[1] %*% t(obs_svd$v[,1,drop=F])
res_mat <- obs_mat - pred_mean_mat

plot(as.numeric(pred_mean_mat), as.numeric(res_mat))
