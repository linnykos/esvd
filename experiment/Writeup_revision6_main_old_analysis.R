rm(list=ls())
load("../results/submission_round1/old_results/step5_clustering_spca.RData")
dim(dat_impute)
svd_res <- svd(dat_impute)
svd_u <- svd_res$u[,1:p] %*% diag(sqrt(svd_res$d[1:p]))


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
opacity_val <- 0.75
col_vec3_svd <- color_func(opacity_val)[c(5, rep(3,2), rep(1,6), rep(2,2),  rep(5,2))]



cluster_center_svd <- .compute_cluster_center(svd_u, .construct_cluster_matrix(cluster_labels))

png(paste0("../../esvd_results/figure/experiment/Revision_writeup6_embedding_comparison.png"),
    height = 1200, width = 2100, res = 300, units = "px")
par(mar = c(5,4,4,2), mfrow = c(1,2))
set.seed(10)
idx <- sample(1:nrow(svd_u))
plot(svd_u[idx,2], -svd_u[idx,3], asp = T, pch = 16, col = col_vec3_svd[cluster_labels][idx],
     cex = 0.75, xlab = "Latent dimension 2", ylab = "Latent dimension 3",
     main = "SVD embedding\n(Constant-variance Gaussian)")

for(ll in 1:nrow(cluster_center)){
  points(cluster_center_svd[ll,i], -cluster_center_svd[ll,j], pch = 16, cex = 2, col = "black")
  points(cluster_center_svd[ll,i], -cluster_center_svd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[[ll]])
}

####

cluster_center_esvd <- .compute_cluster_center(res_our$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))

num_order_vec_esvd <- c(5, rep(3,2), 3, rep(1,3), rep(4,2), rep(2,2),  rep(5,2))
col_vec_esvd <- color_func(1)[num_order_vec_esvd]
col_vec3_esvd <- color_func(opacity_val)[num_order_vec_esvd]

i <- 2; j <- 3
plot(x = res_our$u_mat[,i], y = res_our$u_mat[,j],
     asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
     col = col_vec3_esvd[cluster_labels], pch = 16,
     main = "eSVD embedding\n(Curved Gaussian)")


for(ll in 1:nrow(cluster_center)){
  points(cluster_center[ll,i], cluster_center[ll,j], pch = 16, cex = 2, col = "black")
  points(cluster_center[ll,i], cluster_center[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[[ll]])
}


graphics.off()


############################
col_vec_short <- color_func(0.9)[c(6,4)]

png(paste0("../../esvd_results/figure/experiment/Revision_writeup6_embedding_comparison_trajectory.png"),
    height = 1200, width = 2100, res = 300, units = "px")
par(mar = c(5,4,4,2), mfrow = c(1,2))
set.seed(10)
i <- 2; j <- 3
idx <- sample(1:nrow(svd_u))
plot(svd_u[idx,2], -svd_u[idx,3], asp = T, pch = 16, col = col_vec3_svd[cluster_labels][idx],
     cex = 0.75, xlab = "Latent dimension 2", ylab = "Latent dimension 3",
     main = "SVD embedding\n(Constant-variance Gaussian)")

for(ll in 1:nrow(cluster_center)){
  points(cluster_center_svd[ll,i], -cluster_center_svd[ll,j], pch = 16, cex = 2, col = "black")
  points(cluster_center_svd[ll,i], -cluster_center_svd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[[ll]])
}

curves <- naive_curves$curves
for(ll in 1:length(curves)){
  ord <- curves[[ll]]$ord
  lines(x = curves[[ll]]$s[ord, i], y = -curves[[ll]]$s[ord, j], col = "white", lwd = 8)
  lines(x = curves[[ll]]$s[ord, i], y = -curves[[ll]]$s[ord, j], col = "black", lwd = 4)
}

####

cluster_center_esvd <- .compute_cluster_center(res_our$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))

num_order_vec_esvd <- c(5, rep(3,2), 3, rep(1,3), rep(4,2), rep(2,2),  rep(5,2))
col_vec_esvd <- color_func(1)[num_order_vec_esvd]
col_vec3_esvd <- color_func(opacity_val)[num_order_vec_esvd]

i <- 2; j <- 3
plot(x = res_our$u_mat[,i], y = res_our$u_mat[,j],
     asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
     col = col_vec3_esvd[cluster_labels], pch = 16,
     main = "eSVD embedding\n(Curved Gaussian)")


for(ll in 1:nrow(cluster_center)){
  points(cluster_center[ll,i], cluster_center[ll,j], pch = 16, cex = 2, col = "black")
  points(cluster_center[ll,i], cluster_center[ll,j], pch = 16, cex = 1.5, col = col_vec_esvd[[ll]])
}

curves <- our_curves$curves
for(ll in 1:length(curves)){
  ord <- curves[[ll]]$ord
  lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 8)
  lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = col_vec_short[ll], lwd = 5)
  lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black",
        lty = 3, lwd = 2)
}
graphics.off()
