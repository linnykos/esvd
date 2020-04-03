rm(list=ls())
load("../results/step5_clustering.RData")

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(210/255, 198/255, 36/255, alpha)) #darker yellow
}
cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)

col_vec_svd <- color_func(1)[c(5, rep(3,2), rep(1,6), rep(2,2),  rep(5,2))]
opacity_val <- 0.5
col_vec3_svd <- color_func(opacity_val)[c(5, rep(3,2), rep(1,6), rep(2,2),  rep(5,2))]
col_name <- c("orange", rep("bluish green", 3), rep("yellow", 3), rep("yellow", 2), rep("skyblue", 2), rep("orange", 2))
custom_cluster_group_list <- list(13, 12, 1, c(10,11), c(2,3), c(4:9))
col_info <- data.frame(name = levels(cell_type_vec),
                       idx = sort(unique(cluster_labels)),
                       level = sapply(1:13, function(y){which(sapply(custom_cluster_group_list, function(x){y %in% x}))}),
                       col_name = col_name,
                       col_code = col_vec_svd)
col_info

cluster_center_svd <- .compute_cluster_center(svd_u, .construct_cluster_matrix(cluster_labels))

dim1 <- 1; dim2 <- 2
plot(x = svd_embedding[,dim1], y = svd_embedding[,dim2],
     asp = T, xlab = paste0("Latent dimension ", dim1), ylab = paste0("Latent dimension ", dim2),
     col = col_vec3_svd[cluster_labels], pch = 16,
     main = ifelse(k == 2, "SVD embedding and trajectories\n(Curved Gaussian)","")
)

x <- svd_embedding[,dim1]; y <- svd_embedding[,dim2]
range_vec <- range(c(x,y))
kde_est <- MASS::kde2d(svd_embedding[,dim1], svd_embedding[,dim2], n = 500, lims = c(range_vec, range_vec), h = 2)
plot(NA,  asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
     main = "SVD embedding\n(Constant-variance Gaussian)",
     xlim = range(x), ylim = range(y))
# image(kde_est, col = grDevices::heat.colors(100, alpha = 0.8), xaxt = "n",
#       yaxt = "n", axes = F, ylab = "", xlab = "")
# points(svd_embedding[,i], svd_embedding[,j], col = rgb(0, 0,0, 0.2), pch = 16)
quant_level <- quantile(kde_est$z, probs = 0.95)
contour(kde_est, add = T, drawlabels = F, lwd = 1, levels = quant_level)
points(x = x, y = y, col = col_vec3_svd[cluster_labels], pch = 16)

par(mfrow = c(2,3))
for(i in 1:6){
  clust_idx <- which(col_info[,"level"] %in% i)
  idx <- which(cluster_labels %in% clust_idx)
  plot(x[idx], y[idx], col = col_info[clust_idx[1], "col_code"], xlim = range(x), ylim = range(y),
       pch = 16, asp = T)
}

png(paste0("../../esvd_results/figure/experiment/Revision_writeup6_svd_table.png"),
    height = 1300, width = 1500, res = 300, units = "px")
plot(NA,  asp = T, xlab = paste0("Latent dimension ", dim1), ylab = paste0("Latent dimension ", dim2),
     main = "SVD embedding\n(Constant-variance Gaussian)",
     xlim = c(-10,10), ylim = range(-y), axes = F)
points(x = x, y = -y, col = col_vec3_svd[cluster_labels], pch = 16)
axis(1); axis(2)

for(i in 3:6){
  if(i == 3){
    clust_idx <- which(col_info[,"level"] %in% c(1:3))
  } else {
    clust_idx <- which(col_info[,"level"] %in% i)
  }
  idx <- which(cluster_labels %in% clust_idx)
  kde_est <- MASS::kde2d(svd_embedding[idx,dim1], -svd_embedding[idx,dim2], n = 500, lims = c(range(x), range(-y)))
  col_val <- col_info[clust_idx[1], "col_code"]
  quant_level <- quantile(kde_est$z, probs = 0.925)
  contour(kde_est, add = T, drawlabels = F, lwd = 5, levels = quant_level, col = "white")
  #contour(kde_est, add = T, drawlabels = F, lwd = 2.5, levels = quant_level, col = "black")
  contour(kde_est, add = T, drawlabels = F, lwd = 3, levels = quant_level, col = col_val)
  contour(kde_est, add = T, drawlabels = F, lwd = 1, levels = quant_level,  lty = 3, col = "black")
}
graphics.off()

#####

i <- 1; j <- 3
plot(x = esvd_embedding$u_mat[,i], y = esvd_embedding$u_mat[,j],
     asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
     col = col_vec3_svd[cluster_labels], pch = 16,
     main = ifelse(k == 2, "eSVD embedding and trajectories\n(Curved Gaussian)","")
)
