rm(list=ls())
load("../results/step1_zeisel_gaussian_fitting.RData")
tmp1 <- svd_embedding
tmp2 <- svd_missing_list
load("../results/step3_zeisel_factorization.RData")
svd_embedding <- tmp1
svd_missing_list <- tmp2

# u_mat <- svd_embedding
# u_mat <- esvd_missing_list[[11]][[1]]$u_mat
u_mat <- esvd_embedding$u_mat
rgl::plot3d(u_mat[,1], u_mat[,2], u_mat[,3], asp = T, col = as.numeric(label_vec),
            main = "Zeisel")

#################

mat <- esvd_embedding$u_mat
cluster_labels <- as.numeric(label_vec)
neigh_size <- determine_minimium_neighborhood_size(mat)
print(neigh_size)
neigh_size_vec <- c(3,4,5)
res_list <- lapply(neigh_size_vec, function(x){
  print(x)
  set.seed(10)
  compute_purity(mat, cluster_labels, neighborhood_size = x, num_samples = 5000)
})

for(x in res_list){
  print(quantile(unlist(x$value_list), probs = seq(0, 0.2, length.out = 11)))
}


mat <- svd_embedding
cluster_labels <- as.numeric(label_vec)
neigh_size <- determine_minimium_neighborhood_size(mat)
print(neigh_size)
neigh_size_vec <- c(3,4,5)
res_list2 <- lapply(neigh_size_vec, function(x){
  print(x)
  set.seed(10)
  compute_purity(mat, cluster_labels, neighborhood_size = x, num_samples = 5000)
})

for(x in res_list2){
  print(quantile(unlist(x$value_list), probs = seq(0, 0.25, length.out = 11)))
}

#######################################

color_palatte <- c(rgb(245, 234, 204, maxColorValue = 255), #yellow
                   rgb(189, 57, 60, maxColorValue = 255)) #red

grDevices::png(filename = paste0("../../esvd_results/figure/main/zeisel_purity.png"),
               height = 850, width = 2300, res = 300,
               units = "px")
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
for(i in 1:3){
  plot(NA, xlim = c(0, 0.25), ylim = c(0,1),
       main = paste0("Purity with neighborhood\nsize of ", neigh_size_vec[i]),
       ylab = "Purity score", xlab = "Quantile")

  # draw grid
  for(x in seq(0, 0.25, by = 0.05)){
    lines(rep(x,2), c(-1e4,1e4), col = "gray", lwd = 0.5, lty = 2)
  }
  for(y in seq(0, 1, by = 0.1)){
    lines(c(-1e4,1e4), rep(y,2), col = "gray", lwd = 0.5, lty = 2)
  }


  x_vec <- seq(0, 0.25, length.out = 11)
  y_vec1 <- quantile(unlist(res_list[[i]]$value_list), probs = x_vec)
  points(x_vec, y_vec1, pch = 16, col = "black", cex = 1.7)
  points(x_vec, y_vec1, pch = 16, col = color_palatte[2], cex = 1.5)
  y_vec2 <- quantile(unlist(res_list2[[i]]$value_list), probs = x_vec)
  points(x_vec, y_vec2, pch = 16, col = "black", cex = 1.2)
  points(x_vec, y_vec2, pch = 16, col = color_palatte[1])


  legend("topleft", c("eSVD embedding", "SVD embedding"), fill=color_palatte[c(2,1)], cex=0.75)
}
graphics.off()


#
# ########
# idx <- c(7, 5, 12, 16, 20, 21)
#
# for(i in idx){
#   print(i)
#   set.seed(10)
#
#   mat <- esvd_missing_list[[i]][[1]]$u_mat
#   cluster_labels <- as.numeric(label_vec)
#   neigh_size <- determine_minimium_neighborhood_size(mat)
#   print(neigh_size)
#   set.seed(10)
#   zz <- compute_purity(mat, cluster_labels, neigh_size)
#   print(zz$avg_val)
# }
#
# for(i in c(1:3)){
#   print(i)
#   set.seed(10)
#
#   mat <- esvd_missing_list[[12]][[i]]$u_mat
#   cluster_labels <- as.numeric(label_vec)
#   neigh_size <- determine_minimium_neighborhood_size(mat)
#   print(neigh_size)
#   set.seed(10)
#   zz <- compute_purity(mat, cluster_labels, neigh_size)
#   print(zz$avg_val)
# }
#
#
#
# mat <- svd_embedding
# cluster_labels <- as.numeric(label_vec)
# neigh_size <- determine_minimium_neighborhood_size(mat)
# print(neigh_size)
# set.seed(10)
# zz <- compute_purity(mat, cluster_labels, neigh_size)
# print(zz$avg_val)
