x1 <- c(-2,-1,0,1,2)
x2 <- c(-1,0,1,-1,1)
x3 <- x1-5*x2
x4 <- -x1

x <- cbind(x1,x2,x3,x4)
x

k <- rankMatrix(x)
res_svd <- svd(x)
u_mat <- res_svd$u[,1:k]%*%diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$v[,1:k]%*%diag(sqrt(res_svd$d[1:k]))


y <- x
y[y <= 0] <- 0
k <- rankMatrix(y)

res_svd <- svd(y)
u_mat <- res_svd$u[,1:k]%*%diag(sqrt(res_svd$d[1:k]))
v_mat <- res_svd$v[,1:k]%*%diag(sqrt(res_svd$d[1:k]))
