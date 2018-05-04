rm(list=ls())
load("../../SOUP/data/zeisel.rda")

dat <- zeisel$counts
dat <- t(apply(dat, 1, function(x){x/sum(x)}))
dat <- 1000*log(dat + 1)

sparsity_vec <- apply(dat, 2, function(x){
  length(which(x!=0))
})

sum_vec <- colSums(dat)

quantile_sum <- quantile(sum_vec, prob = 0.95)

idx <- which(sum_vec <= quantile_sum)
dat <- dat[,idx]

anova_func <- function(x){
  uniq <- unique(zeisel$cell.info[,2])
  k <- length(uniq)
  n <- length(x)

  between <- sum(sapply(uniq, function(k){
    idx <- which(zeisel$cell.info[,2] == k)
    length(idx)*(mean(x[idx]) - mean(x))^2
  }))/(k-1)

  within <- sum(sapply(unique(zeisel$cell.info[,2]), function(k){
    idx <- which(zeisel$cell.info[,2] == k)
    sum((x[idx]-mean(x[idx]))^2)
  }))/(n-k)

  between/within
}

anova_vec <- apply(dat, 2, anova_func)

idx <- order(anova_vec, decreasing = T)[1:700]

dat2 <- dat[,idx]
res_svd <- svd(dat2)

dat <- res_svd$u[,1:8]
mat <- .lasso_connectivity(dat, sparsity = 3)
adj <- .coef_to_adj(mat)
g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")

#figure out which samples have no neighbors
vec <- colSums(adj)
idx <- which(vec == 0)
n <- nrow(res_svd$u)
col_vec <- rep(rgb(0,0,0,0.1), n)
col_vec[idx] <- "red"
plot(res_svd$u[,1], res_svd$u[,2], pch = 16, col = col_vec,
     xlim = range(c(res_svd$u[,1], 0)), ylim = range(c(res_svd$u[,2], 0)))


###########################

eig <- eigen(adj)

num_clust <- 8
eigen_vec <- eig$vectors[,1:num_clust]
clustering <- stats::kmeans(eigen_vec,centers = num_clust)

# ploting
col_vec <- c(1:num_clust)[clustering$cluster]
plot(res_svd$u[,1], res_svd$u[,2], pch = 16, col = col_vec,
     xlim = range(c(res_svd$u[,1], 0)), ylim = range(c(res_svd$u[,2], 0)))

#############################
set.seed(10)
l <- igraph::layout.auto(g)
igraph::plot.igraph(g, vertex.label = NA, main = "", vertex.size = 1, layout = l)
