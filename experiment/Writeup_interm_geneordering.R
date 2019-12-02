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
# order_cell <- unlist(pointer_list)

# find the ordering of the cells
tmp_mat <- matrix(NA, nrow = length(which(cluster_labels %in% intersect_labels)), ncol = 2)
tmp_mat[,1] <- which(cluster_labels %in% intersect_labels)
for(i in 1:nrow(tmp_mat)){
  tmp <- sapply(1:nrow(curve), function(x){.l2norm(res_our$u_mat[tmp_mat[i,1], 1:3] - curve[x,])})
  tmp_mat[i,2] <- which.min(tmp)
}
order_cell <- tmp_mat[order(tmp_mat[,2]),1]

# reform the data
tmp_mat <- log(dat_impute[order_cell,]+1)

# apply the changepoint detection
ecp_res <- ecp::e.divisive(X = tmp_mat, sig.lvl=0.05, R=199, k=NULL, min.size=30, alpha=1)

# order the genes just for visualization purposes
hclust_res <- hclust(dist(t(tmp_mat)))
order_gene <- hclust_res$order

clockwise90 = function(a) { t(a[nrow(a):1,]) }
tmp_mat <- log(dat_impute[order_cell, order_gene]+1)
bool_mat <- tmp_mat > 1.5

#png("../figure/experiment/interm_ordering.png", height = 1500, width = 4900, res = 300, units = "px")
png("../figure/experiment/interm_ordering.png", height = 1500, width = 1325, res = 300, units = "px")
image(clockwise90(t(bool_mat)),
      #asp = length(order_gene)/length(order_cell),
      asp = 1,
      col = c("yellow",rgb(55,44,161,maxColorValue = 255)),
      breaks = c(-.5, 0.5, 1.5), xlab = "Cells in pseudo-time ordering", ylab = "Genes",
      xaxt = "n", yaxt = "n")
axis(1, at= seq(100,1000, length.out = 10)/nrow(bool_mat), labels= seq(100,1000, length.out = 10))
axis(2, at= seq(20,200, length.out = 10)/ncol(bool_mat), labels= seq(20,200, length.out = 10))

# draw in all the changepoints now
for(i in 2:(length(ecp_res$estimates)-1)){
  lines(rep(ecp_res$estimates[i]/nrow(bool_mat), 2), c(0,1), lwd = 3, col = "white")
  lines(rep(ecp_res$estimates[i]/nrow(bool_mat), 2), c(0,1), lwd = 2, lty = 2)
}

graphics.off()
