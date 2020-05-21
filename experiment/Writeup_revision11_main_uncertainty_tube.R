rm(list=ls())
load("../results/step5_trajectory_original.RData")

zz <- esvd_sd_res$mat_list
quantile(apply(zz[[1]], 1, function(x){quantile(x, probs = 0.95)}), probs = 0.95)
quantile(apply(zz[[2]], 1, function(x){quantile(x, probs = 0.95)}), probs = 0.95)

x <- apply(zz[[1]], 1, function(x){quantile(x, probs = 0.95)})
plot(x)
x <- apply(zz[[2]], 1, function(x){quantile(x, probs = 0.95)})
plot(x)

zz <- svd_sd_val$mat_list
quantile(apply(zz[[1]], 1, function(x){quantile(x, probs = 0.95)})[200:900], probs = 0.95)
quantile(apply(zz[[2]], 1, function(x){quantile(x, probs = 0.95)})[200:900], probs = 0.95)

x <- apply(zz[[1]], 1, function(x){quantile(x, probs = 0.95)})
plot(x)
x <- apply(zz[[2]], 1, function(x){quantile(x, probs = 0.95)})
plot(x)

#################################
#################################
#################################

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(100/255, 100/255, 100/255, alpha)) #gray
}
color_name_vec <- c("yellow", "skyblue", "bluish green", "blue", "orange", "gray")

#####################

# info for the upcoming svd plots

num_order_vec_svd <- c(5, rep(3,2), c(1,1,1,1,1,1), rep(2,2),  rep(5,2))
col_vec_svd <- color_func(1)[num_order_vec_svd]
col_vec2_svd <- color_func(0.5)[num_order_vec_svd]
col_name_svd <- color_name_vec[num_order_vec_svd]
order_vec_svd <- c(3, 5.1, 5.2, seq(6.1, 6.6, by = 0.1), 4.1, 4.2, 2, 1)
col_info_svd <- data.frame(name = levels(cell_type_vec),
                           idx = sort(unique(cluster_labels)),
                           order = order_vec_svd,
                           col_name = col_name_svd,
                           col_code = col_vec_svd)
col_info_svd$factor_idx <- as.numeric(as.factor(col_info_svd$col_name))
col_info_svd[,c(5,6)] <- col_info_svd[,c(6,5)]
colnames(col_info_svd)[c(5,6)] <- colnames(col_info_svd)[c(6,5)]
col_info_svd
plotting_order_svd <- c(3,1,2,4)

cluster_center_svd <- eSVD:::.compute_cluster_center(svd_embedding[,1:3], .construct_cluster_matrix(cluster_labels))
combn_mat <- utils::combn(3,2)
##################
kk <- 2
graphics::par(mfrow = c(1,3), mar = c(4,4,4,1))

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  graphics::plot(NA, xlim = range(svd_embedding[,i]), ylim = range(svd_embedding[,j]),
                 asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
                 main = ifelse(k == 2, "SVD embedding and trajectories\n(Constant-variance Gaussian)", ""))

  for(ll in plotting_order_svd) {
    target_indices <- col_info_svd$idx[which(col_info_svd$factor_idx %in% ll)]
    idx <- which(cluster_labels %in% target_indices)
    graphics::points(x = svd_embedding[idx,i], y = svd_embedding[idx,j], pch = 16,
                     col = col_vec2_svd[cluster_labels[idx]])
  }

  curves <- svd_bootstrap_list[[kk]]
  for(ll in 1:length(curves)) {
    ord <- curves[[ll]]$ord
    graphics::lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "white", lwd = 5)
    graphics::lines(x = curves[[ll]]$s[ord, i], y = curves[[ll]]$s[ord, j], col = "black", lwd = 3)
  }

  for(ll in 1:nrow(cluster_center_svd)){
    graphics::points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 2, col = "black")
    graphics::points(cluster_center_svd[ll,i], cluster_center_svd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }
}
