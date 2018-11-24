rm(list=ls())
load("../simulation/fit.RData")

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #blue
             rgb(100,100,200,maxColorValue=255), #purple
             rgb(149,219,144,maxColorValue=255)) #green

plot(fit$u_mat[,1], fit$u_mat[,2], asp = T, pch = 16,
     col = col_vec[rep(1:4, each = 50)])

# let's try to fix this
.get_lineages(fit$u_mat, cluster_labels = rep(1:4, each = 50), starting_cluster = 1)

###################

dat <- fit$u_mat
cluster_labels = rep(1:4, each = 50)
starting_cluster = 1
knn <- NA

cluster_mat <- .construct_cluster_matrix(cluster_labels)
k <- ncol(cluster_mat)
stopifnot(is.matrix(dat), nrow(dat) == nrow(cluster_mat))

### get the connectivity matrix
centers <- .compute_cluster_center(dat, cluster_mat)
dat_augment <- rbind(centers, dat)

### construct the k-nearest neighbor graph
if(is.na(knn)){
  knn <- 1
  while(TRUE){
    knn_graph <- .construct_knn_graph(dat_augment, knn = knn)
    if(igraph::components(knn_graph)$no == 1) break()
    knn <- knn + 1
  }
} else {
  knn_graph <- .construct_knn_graph(dat_augment, knn = knn)
  stopifnot(igraph::components(knn_graph)$no == 1)
}

zz <- igraph::betweenness(knn_graph)

# remove outliers based on IQR
iqr_mat <- do.call(rbind, lapply(1:4, function(x){
  idx <- which(cluster_labels == x)
  iqr_val <- stats::IQR(zz[idx+4]); med_val <- stats::median(zz[idx+4])
  cbind(idx, x, (zz[idx+4]-med_val)/iqr_val)
}))
iqr_mat <- iqr_mat[order(iqr_mat[,ncol(iqr_mat)], decreasing = T),]

idx <- iqr_mat[-c(1:floor(nrow(iqr_mat)/20)),1] #remove top 5%

plot(fit$u_mat[idx,1], fit$u_mat[idx,2], asp = T, pch = 16,
     col = col_vec[rep(1:4, each = 50)][idx])
plot(fit$u_mat[,1], fit$u_mat[,2], asp = T, pch = 16,
     col = col_vec[rep(1:4, each = 50)])

# should only remove the top few within each cluster
# and then do the bootstrapping?
# or only keep the edges where the nearest goes both ways?


