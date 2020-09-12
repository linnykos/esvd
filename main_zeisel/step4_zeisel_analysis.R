rm(list=ls())
load("../results/step1_zeisel_gaussian_fitting.RData")
tmp1 <- svd_embedding
tmp2 <- svd_missing_list
load("../results/step3_zeisel_factorization.RData")
svd_embedding <- tmp1
svd_missing_list <- tmp2

#################

mat <- esvd_embedding$u_mat
cluster_labels <- as.numeric(label_vec)
neigh_size <- determine_minimium_neighborhood_size(mat)
neigh_size_vec <- c(3,4,5)
res_list <- lapply(neigh_size_vec, function(x){
  set.seed(10)
  compute_purity(mat, cluster_labels, neighborhood_size = x, num_samples = 5000)
})


vec <- unlist(res_list[[1]]$value_list)
mean(vec)
mean(sort(vec)[1:round(length(vec)/4)])


mat <- svd_embedding
cluster_labels <- as.numeric(label_vec)
neigh_size <- determine_minimium_neighborhood_size(mat)
neigh_size_vec <- c(3,4,5)
res_list2 <- lapply(neigh_size_vec, function(x){
  set.seed(10)
  compute_purity(mat, cluster_labels, neighborhood_size = x, num_samples = 5000)
})

vec <- unlist(res_list2[[1]]$value_list)
mean(vec)
mean(sort(vec)[1:round(length(vec)/4)])

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
