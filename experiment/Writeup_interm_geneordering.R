rm(list=ls())
load("../results/old_results/step5_clustering_spca.RData")

# extract the curve in question
curve <- our_curves$curves$Curve1$s
nrow(dat_impute) == length(cluster_labels)
intersect_labels <- intersect(our_curves$lineages[[1]], our_curves$lineages[[2]])
# pointer_list <- lapply(intersect_labels, function(x){
#   print(paste0("Working on index ", x))
#   idx <- which(cluster_labels == x)
#
#   tmp_mat <- matrix(NA, nrow = length(idx), ncol = 3)
#   tmp_mat[,1] <- idx
#
#   # for each point, find the closest index on the curve
#   for(i in 1:length(idx)){
#     tmp <- sapply(1:nrow(curve), function(x){.l2norm(res_our$u_mat[idx[i], 1:3] - curve[x,])})
#     tmp_mat[i,2] <- which.min(tmp)
#     tmp_mat[i,3] <- min(tmp)
#   }
#
#   tmp_mat[order(tmp_mat[,2]),1]
# })

tmp_mat <- matrix(NA, nrow = length(which(cluster_labels %in% intersect_labels)), ncol = 2)
tmp_mat[,1] <- which(cluster_labels %in% intersect_labels)
for(i in 1:nrow(tmp_mat)){
  tmp <- sapply(1:nrow(curve), function(x){.l2norm(res_our$u_mat[tmp_mat[i,1], 1:3] - curve[x,])})
  tmp_mat[i,2] <- which.min(tmp)
}

order_cell <- tmp_mat[order(tmp_mat[,2]),1]
hclust_res <- hclust(dist(t(dat_impute)))
order_gene <- hclust_res$order

clockwise90 = function(a) { t(a[nrow(a):1,]) }
tmp_mat <- t(log(dat_impute[order_cell, order_gene]+1))
bool_mat <- tmp_mat > 1
image(clockwise90(bool_mat), asp = length(order_gene)/length(order_cell))
