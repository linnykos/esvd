load("../results/step3_factorization_logged.RData")

u_mat <- res$u_mat[,1:4]
dist_vec <- as.numeric(dist(u_mat))
cutoff <- quantile(dist_vec, probs = 0.0025)

clustering <- dbscan::dbscan(u_mat, eps = cutoff)
col_vec <- rep(rgb(0,0,0,0.1), nrow(u_mat))
for(i in 1:max(clustering$cluster)){
  col_vec[which(clustering$cluster == i)] <- i
}

plot(u_mat[,1], u_mat[,2], pch = 16, col = col_vec, asp = T)
plot(u_mat[,1], u_mat[,3], pch = 16, col = col_vec, asp = T)
plot(u_mat[,2], u_mat[,3], pch = 16, col = col_vec, asp = T)
plot(u_mat[,3], u_mat[,4], pch = 16, col = col_vec, asp = T)
