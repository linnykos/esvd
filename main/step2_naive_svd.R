zz <- svd(dat)
k <- 4
u_mat <- zz$u[,1:k] %*% diag(sqrt(zz$d[1:k]))
plot(u_mat[,1], u_mat[,2], pch = 16)
