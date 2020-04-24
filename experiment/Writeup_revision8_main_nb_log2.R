rm(list=ls())
load("../results/step5_trajectory_nb_hvg_log2.RData")

esvd_angle_res
cbind(esvd_angle_res$all_results, paramMat_esvd)

#########

color_func <- function(alpha = 0.2){
  c(rgb(240/255, 228/255, 66/255, alpha), #yellow
    rgb(86/255, 180/255, 233/255, alpha), #skyblue
    rgb(0/255, 158/255, 115/255, alpha), #bluish green
    rgb(0/255, 114/255, 178/255,alpha), #blue
    rgb(230/255, 159/255, 0/255,alpha), #orange
    rgb(100/255, 100/255, 100/255, alpha)) #gray
}

cell_type_vec <- as.character(marques$cell.info$cell.type)
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

###

num_order_vec_svd <- c(5, rep(3,2), c(1,1,1,1,1,1), rep(2,2),  rep(5,2))
col_vec_svd <- color_func(1)[num_order_vec_svd]
col_vec2_svd <- color_func(0.5)[num_order_vec_svd]
col_name_svd <- c("orange", rep("bluish green", 2), rep("yellow", 6), rep("skyblue", 2), rep("orange", 2))
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
plotting_order_svd <- c(2,3,1,4)

cluster_center_svd <- .compute_cluster_center(svd_embedding[,1:3], .construct_cluster_matrix(cluster_labels))
cluster_center_esvd <- .compute_cluster_center(esvd_embedding$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))
combn_mat <- combn(3,2)

################

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision8_nb_2dplots_", k, ".png"),
      height = 1500, width = 1500, res = 300,
      units = "px")
  plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding and trajectories\n(Negative binomial)")

  for(ll in plotting_order_svd) {
    target_indices <- col_info_svd$idx[which(col_info_svd$factor_idx == ll)]
    idx <- which(cluster_labels %in% target_indices)
    points(x = esvd_embedding$u_mat[idx,i], y = esvd_embedding$u_mat[idx,j], pch = 16,
           col = col_info_svd$col_code[target_indices[1]])
  }

  for(ll in 1:nrow(cluster_center_esvd)){
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }
  graphics.off()
}

##################################

# plot each cell type individually (so 13 different plots)
for(zz in unlist(cluster_group_list)){
  png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision8_nb_2dplots_subtype", zz, ".png"),
      height = 800, width = 2500, res = 300,
      units = "px")
  par(mfrow = c(1,3))

  for(k in 1:3){
    i <- combn_mat[1,k]; j <- combn_mat[2,k]

    plot(NA, xlim = range(esvd_embedding$u_mat[,i]), ylim = range(esvd_embedding$u_mat[,j]),
         asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
         main = paste0("For cell type ", col_info_svd[zz, "name"]))

    for(ll in plotting_order_svd) {
      idx <- which(cluster_labels == zz)
      points(x = esvd_embedding$u_mat[idx,i], y = esvd_embedding$u_mat[idx,j], pch = 16,
             col = col_info_svd$col_code[zz])
    }

    for(ll in 1:nrow(cluster_center_esvd)){
      points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
      points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
    }

    points(cluster_center_esvd[zz,i], cluster_center_esvd[zz,j], pch = 16, cex = 2.25*2, col = "black")
    points(cluster_center_esvd[zz,i], cluster_center_esvd[zz,j], pch = 16, cex = 1.5*2, col = col_vec_svd[zz])
  }
  graphics.off()
}

########################################

png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision8_nb_training_testing.png"),
    height = 1500, width = 2500, res = 300,
    units = "px")
par(mfrow = c(1,2))

nat_mat_list <- lapply(1:length(esvd_missing_list[[esvd_angle_res$idx]]), function(i){
  esvd_missing_list[[esvd_angle_res$idx]][[i]]$u_mat %*% t(esvd_missing_list[[esvd_angle_res$idx]][[i]]$v_mat)
})

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 family = "neg_binom",
                                 scalar = paramMat_esvd[esvd_angle_res$idx, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 max_points = 1e6)


plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 family = "neg_binom",
                                 scalar = paramMat_esvd[esvd_angle_res$idx, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Testing set)")

graphics.off()
