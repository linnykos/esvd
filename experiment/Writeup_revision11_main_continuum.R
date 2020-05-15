rm(list=ls())
load("../results/step5_trajectory_original.RData")

session_info2 <- sessionInfo()
source_code_info2 <- ""
date_of_run2 <- Sys.time()

p <- 3
set.seed(10)
esvd_curves <- eSVD::slingshot(esvd_embedding$u_mat[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                               cluster_group_list = cluster_group_list, shrink = 2,
                               verbose = T, upscale_factor = 1, stretch = 9999, max_iter = 3,
                               squared = T)

##########

pseudotime_df <- eSVD:::construct_pseudotime_trajectory_matrix(esvd_curves, cluster_labels)

pseudotime_df2 <- pseudotime_df[-intersect(which(is.na(pseudotime_df$consensus)), which(pseudotime_df$pseudotime <= 3)),]
pseudotime_df2 <- pseudotime_df2[-which(!pseudotime_df2$consensus),]
pseudotime_df2 <- pseudotime_df2[-intersect(which(pseudotime_df2$dist_to_curve >= 0.1),
                                            which(pseudotime_df2$status == "1")),]
traj1_cluster <- c(4,7,6,5)
traj2_cluster <- c(9)

tmp <- pseudotime_df2[which(!pseudotime_df2$cluster_labels %in% traj2_cluster),]
order_vec <- order(tmp$pseudotime, decreasing = F)
x1 <- tmp$pseudotime[order_vec]
dat_ordered1 <- dat_impute[tmp$cell_idx[order_vec],]
idx_trajectory1 <- which(cluster_labels[tmp$cell_idx[order_vec]] %in% traj1_cluster)

tmp <- pseudotime_df2[which(!pseudotime_df2$cluster_labels %in% traj1_cluster),]
order_vec <- order(tmp$pseudotime, decreasing = F)
x2 <- tmp$pseudotime[order_vec]
dat_ordered2 <- dat_impute[tmp$cell_idx[order_vec],]
idx_trajectory2 <- which(cluster_labels[tmp$cell_idx[order_vec]] %in% traj2_cluster)

max_common_idx <- min(c(idx_trajectory1, idx_trajectory2))-1

##############

ncores <- 20
doMC::registerDoMC(cores = ncores)

func <- function(j){
  print(j)
  eSVD:::.find_highly_expressed_region(common_vec = dat_ordered1[1:max_common_idx,j],
                                       specific_vec1 = dat_ordered1[idx_trajectory1,j],
                                       specific_vec2 = dat_ordered1[idx_trajectory2,j],
                                       standardize = T)
}

segmentation_res <- foreach::"%dopar%"(foreach::foreach(j = 1:ncol(dat_impute)), func(j))
save.image("../results/Writeup_revision11_continuum.RData")

