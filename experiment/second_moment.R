mean_vec <- rep(10000,10)
dat <- MASS::mvrnorm(10000, mean_vec, diag(10))
est_mean_vec <- apply(dat, 2, mean)

zz <- mean_vec %*% t(mean_vec) - est_mean_vec %*% t(est_mean_vec)
eigen(zz)
plot(eigen(zz)$values)

.l2norm <- function(x){sqrt(sum(x^2))}
max_val <- max(abs(eigen(zz)$values))
.l2norm(mean_vec - est_mean_vec) * .l2norm(mean_vec + est_mean_vec)

max(abs(eigen(zz)$values))
.l2norm(mean_vec - est_mean_vec)
sqrt(t(mean_vec)%*%est_mean_vec)

##############

set.seed(10)
mean_vec <- rep(10000,10)
n <- 10000
dat <- MASS::mvrnorm(n, mean_vec, diag(10))
est_mean_vec <- apply(dat, 2, mean)

est_moment <- t(dat) %*% dat/n
moment <- mean_vec %*% t(mean_vec)
zz <- moment - est_moment
plot(eigen(zz)$values)
max(eigen(zz)$values)

est_cov <- stats::cov(dat)
cov <- diag(10)
zz_diff <- est_cov - cov
max(eigen(zz_diff)$values)

#############

svd_decomp <- svd(dat)
sum(abs(dat - svd_decomp$u %*% diag(svd_decomp$d) %*% t(svd_decomp$v))) < 1e-6
u <- svd_decomp$u
lambda <- svd_decomp$d^2/n
v <- svd_decomp$v

sum(abs(svd_decomp$u %*% diag(svd_decomp$d) %*% t(svd_decomp$v) -
          sqrt(n)*u %*% diag(sqrt(lambda)) %*% t(v)))

est_moment2 <- v%*%diag(lambda)%*%t(v)
zz <- moment - est_moment2
max(eigen(zz)$values)
