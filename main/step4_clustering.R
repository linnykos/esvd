load("../results/step3_factorization_logged.RData")

u_mat <- res$u_mat[,1:3]

# custom dbscan
dist_mat <- as.matrix(dist(u_mat))
min_dist <- quantile(apply(dist_mat, 1, function(x){
  quantile(x, 0.05)
}), 0.75)
n <- nrow(u_mat)
neighbor_count <- 10
adj_mat <- matrix(0, n, n)
for(i in 1:n){
  idx1 <- which(dist_mat[,i] <= min_dist)
  idx2 <- order(dist_mat[,i], decreasing = F)[1:neighbor_count]
  idx <- intersect(idx1, idx2)
  if(length(idx) > 0) adj_mat[idx,i] <- 1
}
adj_mat <- adj_mat + t(adj_mat)
adj_mat[which(adj_mat > 0)] <- 1

# create the graph
neigh_vec <- apply(adj_mat, 2, function(x){length(which(x != 0))})
g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
g <- igraph::simplify(g)
assign_vec <- rep(NA, n)
upper_cutoff <- 14

count_idx <- 0
for(i in 1:n){
  if(!is.na(assign_vec[i])) next()
  if(neigh_vec[i] < upper_cutoff){assign_vec[i] <- 0; next()}

  q <- dequer::queue()
  count_idx <- count_idx + 1
  assign_vec[i] <- count_idx
  tmp <- sapply(as.numeric(igraph::neighbors(g, i)), function(x){
    dequer::pushback(q, x); invisible()
  })

  while(length(q) > 0){
    item <- dequer::pop(q)
    if(!is.na(assign_vec[item]) & assign_vec[item] != 0) next()
    assign_vec[item] <- count_idx
    if(neigh_vec[item] >= upper_cutoff){
      tmp <- sapply(as.numeric(igraph::neighbors(g, item)), function(x){
        dequer::pushback(q, x); invisible()
      })
    }
  }
}

# remove small clusters
size_cutoff <- 19
tab <- table(assign_vec)
cluster_remove <- as.numeric(names(tab)[which(tab <= size_cutoff)])
assign_vec[assign_vec %in% cluster_remove] <- 0
mat <- cbind(sort(unique(assign_vec))[-1], 1:(length(unique(assign_vec))-1))
for(i in 1:nrow(mat)){
  assign_vec[which(assign_vec == mat[i,1])] <- mat[i,2]
}
table(assign_vec)

col_vec <- assign_vec
col_vec[which(col_vec == 0)] <- rgb(0,0,0,0.1)
plot(u_mat[,1], u_mat[,2], pch = 16, col = col_vec, asp = T)
plot(u_mat[,1], u_mat[,3], pch = 16, col = col_vec, asp = T)
plot(u_mat[,2], u_mat[,3], pch = 16, col = col_vec, asp = T)

####################

# specific clustering
dist_vec <- as.numeric(dist(u_mat))
clustering <- dbscan::dbscan(u_mat, eps = 0.1)
clustering

cutoff <- quantile(dist_vec, probs = 0.002)
clustering <- dbscan::dbscan(u_mat, eps = cutoff)

# general clustering
dist_vec2 <- as.numeric(dist(u_mat[,1:2]))
cutoff2 <- quantile(dist_vec2, probs = 0.025)
clustering2 <- dbscan::dbscan(u_mat[,1:2], eps = cutoff2)

assignment_vec <- rep(NA, nrow(u_mat))
for(i in 1:max(clustering$cluster)){
  idx <- which(clustering$cluster == i)
  assignment_vec[idx] <- i
}
for(i in 2:max(clustering2$cluster)){
  idx <- which(clustering2$cluster == i)
  stopifnot(all(is.na(assignment_vec[idx])))
  assignment_vec[idx] <- i + max(clustering$cluster) - 1
}

# remove small clusters
idx <- which(table(assignment_vec) <= 10)
assignment_vec[which(assignment_vec %in% idx)] <- NA
mat <- cbind(sort(unique(assignment_vec)), 1:c(length(unique(assignment_vec))-1))
for(i in 1:nrow(mat)){
  assignment_vec[which(assignment_vec == mat[i,1])] <- mat[i,2]
}

###########

col_vec <- assignment_vec
col_vec[which(is.na(col_vec))] <- rgb(0,0,0,0.1)

plot(u_mat[,1], u_mat[,2], pch = 16, col = col_vec, asp = T)
plot(u_mat[,1], u_mat[,3], pch = 16, col = col_vec, asp = T)
# plot(u_mat[,1], u_mat[,4], pch = 16, col = col_vec, asp = T)
plot(u_mat[,2], u_mat[,3], pch = 16, col = col_vec, asp = T)
# plot(u_mat[,2], u_mat[,4], pch = 16, col = col_vec, asp = T)
# plot(u_mat[,3], u_mat[,4], pch = 16, col = col_vec, asp = T)

#########

# rearrange dat_impute
dat2 <- dat_impute_log
idx <- which(is.na(assignment_vec))
assign_max <- max(assignment_vec, na.rm = T)
for(i in assign_max:1){
  idx <- c(which(assignment_vec == i), idx)
}

.plot_singlecell(dat2[idx,])
line_idx <- sapply(1:assign_max, function(x){
  1-length(which(assignment_vec <= x))/length(assignment_vec)
})
for(i in 1:assign_max){
  lines(x = c(0,1), y = rep(line_idx[i], 2), lwd = 2, lty = 2)
}

#######

# compare this clustering to kmeans
clustering_kmeans <- kmeans(dat_impute_log, centers = 10, nstart = 10, iter.max = 50)
idx <- as.numeric(unlist(lapply(1:10, function(x){
  which(clustering_kmeans$cluster == x)
})))
.plot_singlecell(dat2[idx,])

# oh that worked really well. let's just paint using this

plot(u_mat[,1], u_mat[,2], pch = 16, col = clustering_kmeans$cluster, asp = T)
plot(u_mat[,1], u_mat[,3], pch = 16, col = clustering_kmeans$cluster, asp = T)
plot(u_mat[,2], u_mat[,3], pch = 16, col = clustering_kmeans$cluster, asp = T)
plot(u_mat[,2], u_mat[,4], pch = 16, col = clustering_kmeans$cluster, asp = T)
plot(u_mat[,3], u_mat[,4], pch = 16, col = clustering_kmeans$cluster, asp = T)

# hm not that good

########################################
#
# set.seed(10)
# zz = MeanShift::bmsClustering(t(u_mat), h = 0.8)
# idx <- which(zz$labels == 1)
# zz2 <- MeanShift::bmsClustering(t(u_mat[idx,]), h = 0.5)
#
# table(zz2$labels)
# plot(u_mat[idx,1], u_mat[idx,2], col = zz2$labels, asp = T, pch = 16)

