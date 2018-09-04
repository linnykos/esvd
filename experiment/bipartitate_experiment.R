rm(list=ls())

# set up parameters
adj <- matrix(c(0, 0.1, -.05,
                0.15,-.1,-.1,
                -.25,-.15,0.05), 3, 3, byrow = T)
tmp <- svd(adj)
u_center <- t(tmp$u %*% diag(sqrt(tmp$d)))
v_center <- t(tmp$v %*% diag(sqrt(tmp$d)))


u_num <- c(20, 10, 80)
u_label <- unlist(lapply(1:3, function(x){rep(x, u_num[x])}))
v_num <- c(30, 40, 100)
v_label <- unlist(lapply(1:3, function(x){rep(x, v_num[x])}))

# generate matrices
set.seed(10)
u_sig <- 0.03
u_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = u_num[x], mu = u_center[,x], Sigma = u_sig*diag(3))
}))
v_sig <- 0.03
v_dat <- do.call(rbind, lapply(1:3, function(x){
  MASS::mvrnorm(n = v_num[x], mu = v_center[,x], Sigma = v_sig*diag(3))
}))

mean_dat <- u_dat %*% t(v_dat)
dim(mean_dat)

all_dat <- rbind(u_dat, v_dat)
all_gram <- all_dat %*% t(all_dat)

# svd
res_svd <- svd(mean_dat)
svd_u <- res_svd$u %*% diag(sqrt(res_svd$d))
svd_v <- res_svd$v %*% diag(sqrt(res_svd$d))
res_svd$d[1:10]

# eigen
res_eig <- eigen(all_gram)
res_eig$values[1:10]

############################

mat = matrix(1:8, ncol = 4, nrow = 2)
res = svd(mat)
mat = res$u %*% diag(2:1) %*% t(res$v)
mat2 = matrix(0, 6, 6)
mat2[1:2, 3:6] = mat
mat2[3:6, 1:2] = t(mat)

svd(mat)
svd(mat2)
eigen(mat2)
svd(mat)

svd(mat %*% t(mat))

