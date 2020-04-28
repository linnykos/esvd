rm(list=ls())
load("../results/step3_scalar_heuristic_cg_hvg_tmp.RData")

##############################

#################

idx_tuned <- 7
png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision8_curved_gaussian_training_testing.png"),
    height = 1500, width = 2500, res = 300,
    units = "px")
par(mfrow = c(1,2))

nat_mat_list <- lapply(1:length(esvd_missing_list[[idx_tuned]]), function(i){
  esvd_missing_list[[idx_tuned]][[i]]$u_mat %*% t(esvd_missing_list[[idx_tuned]][[i]]$v_mat)
})

training_idx_list <- lapply(1:length(missing_idx_list), function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = training_idx_list,
                                 family = "curved_gaussian",
                                 scalar = paramMat_esvd[idx_tuned, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Training set)",
                                 max_points = 1e6)


plot_prediction_against_observed(dat_impute, nat_mat_list = nat_mat_list,
                                 missing_idx_list = missing_idx_list,
                                 family = "curved_gaussian",
                                 scalar = paramMat_esvd[idx_tuned, "scalar"],
                                 main = "eSVD embedding:\nMatrix-completion diagnostic\n(Testing set)")

graphics.off()

###############################3

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
cluster_center_esvd <- .compute_cluster_center(esvd_missing_list[[idx_tuned]][[1]]$u_mat[,1:3], .construct_cluster_matrix(cluster_labels))
combn_mat <- combn(3,2)

################


for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  # png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision8_curved_gaussian_2dplots_", k, ".png"),
  #     height = 1500, width = 1500, res = 300,
  #     units = "px")
  plot(x = esvd_missing_list[[idx_tuned]][[1]]$u_mat[,i],
       y = esvd_missing_list[[idx_tuned]][[1]]$u_mat[,j],
       xlim = range(esvd_missing_list[[idx_tuned]][[1]]$u_mat[,i]), ylim = range(esvd_missing_list[[idx_tuned]][[1]]$u_mat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding\n(Curved gaussian)",
       col = col_info_svd[cluster_labels,"col_code"], pch = 16)

  for(ll in 1:nrow(cluster_center_esvd)){
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }
  # graphics.off()
}

zz <- esvd_missing_list[[idx_tuned]][[1]]$u_mat
tmp <- stats::prcomp(zz, center = T, scale. = F)
tmp_dat <- tmp$x[,1:3]
cluster_center_esvd <- .compute_cluster_center(tmp_dat, .construct_cluster_matrix(cluster_labels))

for(k in 1:ncol(combn_mat)){
  i <- combn_mat[1,k]; j <- combn_mat[2,k]

  # png(filename = paste0("../../esvd_results/figure/experiment/Writeup_revision8_curved_gaussian_2dplots_", k, "_janky.png"),
  #     height = 1500, width = 1500, res = 300,
  #     units = "px")
  plot(x = tmp_dat[,i],
       y = tmp_dat[,j],
       xlim = range(tmp_dat[,i]), ylim = range(tmp_dat[,j]),
       asp = T, xlab = paste0("Latent dimension ", i), ylab = paste0("Latent dimension ", j),
       main = "eSVD embedding\n(Curved gaussian)",
       col = col_info_svd[cluster_labels,"col_code"], pch = 16)

  for(ll in 1:nrow(cluster_center_esvd)){
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 2.25, col = "black")
    points(cluster_center_esvd[ll,i], cluster_center_esvd[ll,j], pch = 16, cex = 1.5, col = col_vec_svd[ll])
  }
  # graphics.off()
}

###############################3

training_idx_list <- lapply(1:cv_trials, function(i){
  c(1:prod(dim(dat_impute)))[-missing_idx_list[[i]]]
})

nat_mat_list_list <- lapply(1:nrow(paramMat_esvd), function(i){
  lapply(1:cv_trials, function(j){
    u_mat <- esvd_missing_list[[i]][[j]]$u_mat
    v_mat <- esvd_missing_list[[i]][[j]]$v_mat
    u_mat %*% t(v_mat)
  })
})

esvd_angle_res <- eSVD:::tuning_select_scalar(dat = dat_impute, nat_mat_list_list = nat_mat_list_list,
                                              family = fitting_distr,  missing_idx_list = training_idx_list,
                                              scalar_vec = paramMat_esvd[,"scalar"])
# okay good, the training results are also pretty ass
