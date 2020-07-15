rm(list=ls())
load("../results/submission_round2/step5_trajectory_original.RData")

res <- prepare_data_for_segmentation(dat_impute, cluster_labels, curve_list = esvd_curves_long,
                                     min_traj_pseudotime = 6, cluster_removal_idx_vec = c(2,3,4),
                                     cluster_removal_time_vec = c(7,7,7))

cell_idx_common <- res$cell_idx_common
cell_idx_traj1 <- res$cell_idx_traj1
cell_idx_traj2 <- res$cell_idx_traj2

name_vec <- ls()
name_vec <- name_vec[-which(name_vec %in% c("dat_impute", "cell_idx_common", "cell_idx_traj1", "cell_idx_traj2"))]
rm(list = name_vec)

save.image("../../../../Jing/dynamicGraphRoot/data-raw/esvd_marques.RData")
