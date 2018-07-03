set.seed(10)
x <- matrix(rnorm(100), 10, 10)
x <- t(x) %*% x
x <- -x

res <- eigen(x)
full_x <- res$vectors %*% diag(res$values) %*% t(res$vectors)

approx_1 <- res$vectors[,9:10] %*% diag(res$values[9:10]) %*% t(res$vectors[,9:10])
sum((approx_1 - x)^2)

approx_2 <- res$vectors[,1:2] %*% diag(res$values[1:2]) %*% t(res$vectors[,1:2])
sum((approx_2 - x)^2)

vec <- res$vectors[,10]
t(vec) %*% x %*% vec

########################

#construct a new matrix
vec <- res$values
vec[9] <- -vec[9]
vec[1:9] <- vec[c(9,1:8)]
mat <- res$vectors
mat[,1:9] <- mat[,c(9,1:8)]
y <- mat %*% diag(vec) %*% t(mat)

res <- eigen(y)
approx_1 <- res$vectors[,c(1,10)] %*% diag(res$values[c(1,10)]) %*% t(res$vectors[,c(1,10)])
sum((approx_1 - y)^2)

approx_2 <- res$vectors[,c(1,2)] %*% diag(res$values[c(1,2)]) %*% t(res$vectors[,c(1,2)])
sum((approx_2 - y)^2)

approx_3 <- res$vectors[,c(9,10)] %*% diag(res$values[c(9,10)]) %*% t(res$vectors[,c(9,10)])
sum((approx_3 - y)^2)

