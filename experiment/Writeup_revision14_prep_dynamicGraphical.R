rm(list=ls())
load("../results/submission_round2/step5_trajectory_original.RData")

res <- prepare_data_for_segmentation(dat_impute, cluster_labels, curve_list = esvd_curves_long,
                                     min_traj_pseudotime = 6, cluster_removal_idx_vec = c(2,3,4),
                                     cluster_removal_time_vec = c(7,7,7))

cell_idx_common <- res$cell_idx_common
cell_idx_traj1 <- res$cell_idx_traj1
cell_idx_traj2 <- res$cell_idx_traj2

########

color_func <- function(alpha = 0.2){
  c(grDevices::rgb(240/255, 228/255, 66/255, alpha), #yellow
    grDevices::rgb(86/255, 180/255, 233/255, alpha), #skyblue
    grDevices::rgb(0/255, 158/255, 115/255, alpha), #bluish green
    grDevices::rgb(0/255, 114/255, 178/255,alpha), #blue
    grDevices::rgb(230/255, 159/255, 0/255,alpha), #orange
    grDevices::rgb(100/255, 100/255, 100/255, alpha)) #gray
}
color_name_vec <- c("yellow", "skyblue", "bluish green", "blue", "orange", "gray")

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

dat_smooth <- 1/(esvd_embedding$u_mat %*% t(esvd_embedding$v_mat))
colnames(dat_smooth) <- colnames(dat_impute)

name_vec <- ls()
name_vec <- name_vec[-which(name_vec %in% c("dat_impute", "dat_smooth", "cell_idx_common", "cell_idx_traj1",
                                            "cell_idx_traj2", "cluster_labels", "col_info_svd"))]
rm(list = name_vec)
rm(list = "name_vec")

save.image("../../../../Jing/dynamicGraphRoot/data-raw/esvd_marques.RData")
