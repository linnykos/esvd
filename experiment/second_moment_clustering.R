d <- 40
n <- 600
k <- 4
tol <- 0.1

c_list <- list(4)
c_list[[1]] <- matrix(c(0.9, 0.9, 0.3, 0,     0.9, 0.9, 0.2, 0.1,
                0.3, 0.2, 0.9, 0.2,           0, 0.1, 0.2, 0.9), 4, 4)
c_list[[2]] <- matrix(c(0.9, 0.2, 0.5, 0.4,   0.2, 0.9, 0.9, 0,
                0.5, 0.9, 0.9, 0.2,           0.4, 0, 0.2, 0.9), 4, 4)
c_list[[3]] <- matrix(c(0.9, 0.5, 0.5, 0.1,   0.5, 0.9, 0.5, 0.3,
                0.5, 0.5, 0.9, 0.9,           0.1, 0.3, 0.9, 0.9), 4, 4)
c_list[[4]] <- matrix(c(0.9, 0, 0, 0,         0, 0.9, 0, 0,
                0, 0, 0.9, 0,                 0, 0, 0, 0.9), 4, 4)

assignment_mat <- do.call("rbind", lapply(1:4, function(x){
  mat <- matrix(0, floor(d/k), 4)
  mat[,x] <- 1
  mat
}))

covar_list <- lapply(c_list, function(x){
  mat <- assignment_mat %*% x %*% t(assignment_mat) + diag(rep(0.1, d))
  eig <- eigen(mat)
  eig$values[eig$values < tol] <- tol
  eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
})

set.seed(10)
tmp <- matrix(rnorm(d^2),d,d)
tmp <- tmp + t(tmp)
basis_vec <- eigen(tmp)$vectors[,1:2]

set.seed(10)
covar_assignment <- rep(c(1:k), times = floor(n/k))
dat <- t(sapply(1:length(covar_assignment), function(x){
  mean_vec <- 5*rnorm(1)*basis_vec[,1] + 5*rnorm(1)*basis_vec[,2]
  MASS::mvrnorm(1, mean_vec, covar_list[[covar_assignment[x]]])
}))

################

svd_res <- svd(dat)
plot(svd_res$d)
col_vec <- c(1:4)[covar_assignment]
plot(svd_res$v[,1], svd_res$v[,2], pch = 16, col = col_vec)
plot(svd_res$u[,1], svd_res$u[,2], pch = 16)

#################

dat_demean <- dat - svd_res$u[,1:2] %*% diag(svd_res$d[1:2]) %*% t(svd_res$v[,1:2])
svd_res2 <- svd(dat_demean)
plot(svd_res2$v[,1], svd_res2$v[,2], pch = 16, col = col_vec)
image(cov(dat_demean[which(covar_assignment == 1),]))
image(cov(dat_demean[which(covar_assignment == 2),]))
image(cov(dat_demean[which(covar_assignment == 3),]))
image(cov(dat_demean[which(covar_assignment == 4),]))

###################

dat_unraveled <- t(apply(dat_demean, 1, function(x){
  mat <- x %*% t(x)
  mat[upper.tri(mat, diag = T)]
}))
svd_res3 <- svd(dat_unraveled)

plot(svd_res3$v[,1], svd_res3$v[,2], pch = 16)
