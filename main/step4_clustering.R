load("../results/step3_factorization_logged.RData")

# specific clustering
u_mat <- res$u_mat[,1:4]
dist_vec <- as.numeric(dist(u_mat))
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
plot(u_mat[,1], u_mat[,4], pch = 16, col = col_vec, asp = T)
plot(u_mat[,2], u_mat[,3], pch = 16, col = col_vec, asp = T)
plot(u_mat[,2], u_mat[,4], pch = 16, col = col_vec, asp = T)
plot(u_mat[,3], u_mat[,4], pch = 16, col = col_vec, asp = T)

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

