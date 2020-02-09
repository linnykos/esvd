rm(list=ls())
set.seed(10)
u_mat <- MASS::mvrnorm(60, rep(0, 5), diag(5))
v_mat <- MASS::mvrnorm(50, rep(1, 5), 2*diag(5))
n <- nrow(u_mat)
p <- nrow(v_mat)

pred_mat <- u_mat %*% t(v_mat)
svd_res <- svd(pred_mat)
u_sol <- (n/p)^(1/4)*svd_res$u[,1:5] %*% diag(sqrt(svd_res$d[1:5]))
v_sol <- (p/n)^(1/4)*svd_res$v[,1:5] %*% diag(sqrt(svd_res$d[1:5]))

t(u_sol) %*% u_sol/n
t(v_sol) %*% v_sol/p

##############################

svd_u <- svd(u_mat)
svd_v <- svd(v_mat)

svd_3 <- svd(.diag_matrix(svd_u$d) %*% t(svd_u$v) %*% svd_v$v %*% .diag_matrix(svd_u$d))

u_atm <- (n/p)^(1/4)*u_mat %*% svd_3$u %*% .diag_matrix(sqrt(svd_3$d))
v_atm <- (p/n)^(1/4)*v_mat %*% svd_3$v %*% .diag_matrix(sqrt(svd_3$d))

zz = u_mat %*% svd_3$u
t(zz) %*% zz

###########################

zz = .identification(t(u_mat)%*%u_mat, t(v_mat)%*%v_mat)
eigen_res <- eigen(zz %*% t(u_mat) %*% u_mat %*% t(zz))


u_atm = t(sapply(1:nrow(u_mat), function(i){
  zz %*% u_mat[i,]
}))
u_atm2 <- u_mat %*% t(zz)
u_atm3 <- u_mat %*% t(zz) %*% eigen_res$vectors

v_atm = t(sapply(1:nrow(v_mat), function(i){
  solve(t(zz)) %*% v_mat[i,]
}))
v_atm2 <- v_mat %*% solve(zz)
v_atm3 <- v_mat %*% solve(zz) %*% eigen_res$vectors

t(u_atm3) %*% u_atm3
t(v_atm3) %*% v_atm3

pred_mat <- u_mat %*% t(v_mat)
pred_mat[1:5,1:5]
pred_mat3 <- u_atm3 %*% t(v_atm3)
pred_mat3[1:5,1:5]
