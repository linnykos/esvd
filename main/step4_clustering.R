load("../results/step3_factorization.RData")

X <- t(res$u_mat[,1:3])
quantile( dist( t( X ) ), 0.3 )
res_cluster <- MeanShift::bmsClustering(t(res$u_mat[,1:3]), h = 0.8)
table(res_cluster$labels)

plot(res$u_mat[,1], res$u_mat[,2], pch = 16, col = c(1:max(res_cluster$labels))[res_cluster$labels], asp = T)


# maybe try DBSCAN or something that can adaptively cluster?

# res_clustering <- SOUP::SOUP(res$u_mat[,1:3], type = "count", Ks = 5, nPC = 5)

# # let's try spectral clustering first
# normalized <- res$u_mat[,1:3]
# normalized <- t(apply(normalized, 1, function(x){x/.l2norm(x)}))
#
# set.seed(10)
# clustering <- stats::kmeans(normalized, centers = 4, iter.max = 50, nstart = 5)
#
# col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
#              rgb(180,200,255,maxColorValue=255), #blue
#              rgb(100,100,200,maxColorValue=255), #purple
#              rgb(149,219,144,maxColorValue=255)) #green
#
# plot(res$u_mat[,1], res$u_mat[,2], pch = 16, col = col_vec[clustering$cluster], asp = T)
