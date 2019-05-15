rm(list=ls())
load("../experiment/tmp.RData")
paramMat <- cbind(50, 120, 0.05, 150, 2, 2, -2000)
colnames(paramMat) <- c("n_each", "d_each", "sigma", "total", "k", "scalar", "max_val")
vec <- paramMat[1,]

plot(res_our[,1], res_our[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = paramMat[1,"n_each"])])

cluster_labels <- rep(1:4, each = vec["n_each"])
res <- slingshot(res_our, cluster_labels, starting_cluster = 1,
                 use_initialization = T, reduction_percentage = 0.2)

plot(res_our[,1], res_our[,2], asp = T, pch = 16, col = c(1:4)[rep(1:4, each = paramMat[1,"n_each"])])
for(i in 1:2){
  ord <- res$curves[[i]]$ord
  lines(res$curves[[i]]$s[ord,1], res$curves[[i]]$s[ord,2], lwd = 2)
}
