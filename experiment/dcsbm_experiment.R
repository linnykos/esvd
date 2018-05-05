load("data/dat.rda")
res_svd <- svd(dat)

plot(res_svd$u[,1], res_svd$u[,2], pch = 16, col = rgb(0.1, 0.1, 0.1, 0.25),
     xlim = range(c(res_svd$u[,1], 0)), ylim = range(c(res_svd$u[,2], 0)))

lines(c(-10,10), rep(0,2), col = "red", lwd = 2)
lines(rep(0,2), c(-10,10), col = "red", lwd = 2)

u_mat <- t(apply(res_svd$u[,1:2], 1, function(x){x/.l2norm(x)}))
plot(u_mat[,1], u_mat[,2], pch = 16, col = rgb(0.1, 0.1, 0.1, 0.1), asp = T)

clust <- kmeans(u_mat, centers = 5, iter.max = 100, nstart = 10)
col_vec <- c(rgb(1,0,0,0.1), rgb(0,1,0,0.1), rgb(0,0,1,0.1), rgb(0,0,0,0.1),
             rgb(1,0,1,0.1))[clust$cluster]
plot(u_mat[,1], u_mat[,2], pch = 16, col = col_vec, asp = T)


plot(res_svd$u[,1], res_svd$u[,2], pch = 16, col = col_vec, xlim = range(c(res_svd$u[,1], 0)), ylim = range(c(res_svd$u[,2], 0)))

load("../../SOUP/data/zeisel.rda")
zz <- table(zeisel$cell.info[,2], clust$cluster)
zz1 <- t(apply(zz, 1, function(x){round(x/sum(x),2)*100}))

###

# try more dimensions
u_mat <- t(apply(res_svd$u[,1:8]%*%diag(res_svd$d[1:8]), 1, function(x){x/.l2norm(x)}))
clust <- kmeans(u_mat, centers = 5, iter.max = 100, nstart = 10)
col_vec <- c(rgb(1,0,0,0.1), rgb(0,1,0,0.1), rgb(0,0,1,0.1), rgb(0,0,0,0.1),
             rgb(1,0,1,0.1))[clust$cluster]
plot(u_mat[,1], u_mat[,2], pch = 16, col = col_vec, asp = T)

plot(res_svd$u[,1]*res_svd$d[1], res_svd$u[,2]*res_svd$d[2], pch = 16, col = col_vec)

zz2 <- table(zeisel$cell.info[,2], clust$cluster)
zz3 <- t(apply(zz2, 1, function(x){round(x/sum(x),2)*100}))
