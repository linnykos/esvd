rm(list=ls())
load("../results/step5_trajectory_original.RData")

session_info2 <- sessionInfo()
source_code_info2 <- ""
date_of_run2 <- Sys.time()

set.seed(10)
esvd_curves <- eSVD::slingshot(esvd_embedding$u_mat, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                               cluster_group_list = cluster_group_list, shrink = 3,
                               verbose = T, upscale_factor = 1, stretch = 9999, max_iter = 3,
                               squared = T)

##########

pseudotime_df <- eSVD:::construct_pseudotime_trajectory_matrix(esvd_curves, cluster_labels)

pseudotime_df2 <- pseudotime_df[-intersect(which(is.na(pseudotime_df$consensus)), which(pseudotime_df$pseudotime <= 6)),]
pseudotime_df2 <- pseudotime_df2[-which(!pseudotime_df2$consensus),]

traj1_cluster <- c(6,5)
traj2_cluster <- c(7,8,9)

pseudotime_df2 <- pseudotime_df2[order(pseudotime_df2$pseudotime),]
tmp1 <- pseudotime_df2[which(!pseudotime_df2$cluster_labels %in% traj2_cluster),]
tmp2 <- pseudotime_df2[which(!pseudotime_df2$cluster_labels %in% traj1_cluster),]
pseudotime_max_common <- min(pseudotime_df2$pseudotime[pseudotime_df2$cluster_labels %in% c(traj1_cluster, traj2_cluster)])
cell_idx_common <- pseudotime_df2$cell_idx[which(pseudotime_df2$pseudotime <= pseudotime_max_common)]

cell_idx_traj1 <- pseudotime_df2$cell_idx[intersect(which(pseudotime_df2$pseudotime >= pseudotime_max_common),
                                                    which(pseudotime_df2$cluster_labels %in% traj1_cluster))]
cell_idx_traj2 <- pseudotime_df2$cell_idx[intersect(which(pseudotime_df2$pseudotime >= pseudotime_max_common),
                                                    which(pseudotime_df2$cluster_labels %in% traj2_cluster))]

dat1 <- dat_impute[c(cell_idx_common,cell_idx_traj1),]
dat2 <- dat_impute[c(cell_idx_common,cell_idx_traj2),]

segmentation_res <- eSVD:::segment_genes_along_trajectories(dat1, dat2, common_n = length(cell_idx_common),
                                                            standardize = T, verbose = T, ncores = 20)

save.image("../results/Writeup_revision11_continuum.RData")

