rm(list=ls())
source("../experiment/Week23_simulation_generator.R")

set.seed(10)
load("../experiment/Week23_tmp.RData")

col_vec <- c(rgb(205,40,54,maxColorValue=255),
             rgb(180,200,255,maxColorValue=255),
             rgb(100,100,200,maxColorValue=255),
             rgb(149,219,144,maxColorValue=255))

length(which(dat == 0))/prod(dim(dat))
.plot_singlecell(dat)

par(mfrow = c(1,3))
zz <- apply(dat, 1, function(x){length(which(x != 0))/length(x)})
plot(sort(zz))
zz <- apply(dat, 2, function(x){length(which(x != 0))/length(x)})
plot(sort(zz))
plot(sort(dat[dat != 0]))


###########

# try the naive thing
svd_res <- svd(dat)
k <- 5
u_mat_naive <- svd_res$u[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
v_mat_naive <- svd_res$v[,1:k] %*% diag(sqrt(svd_res$d[1:k]))
library(slingshot)
slingshot_res <- slingshot::slingshot(data = u_mat_naive, clusterLabels = rep(1:simulation$h, each = simulation$n_each))

plot(u_mat_naive[,1], u_mat_naive[,2],
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Estimated dim. 1", ylab = "Estimated dim. 2", main = "Cell estimated vectors")
lines(slingshot_res)

# let's see how well this naive method did
.evaluate_objective(dat, u_mat_naive, v_mat_naive)

##############

# compare to our implemented method
.evaluate_objective(dat, res_nodropout$u_mat, res_nodropout$v_mat)

plot(res_nodropout$u_mat[,1], res_nodropout$u_mat[,2],
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Estimated dim. 1", ylab = "Estimated dim. 2", main = "Cell estimated vectors")

