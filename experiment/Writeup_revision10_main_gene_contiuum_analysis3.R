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

######################

ncores <- 20
doMC::registerDoMC(cores = ncores)

func <- function(j){
  print(j)
  eSVD:::.find_highly_expressed_region(dat_ordered1[,j], idx_trajectory1,
                                dat_ordered2[,j], idx_trajectory2)
}

segmentation_res <- foreach::"%dopar%"(foreach::foreach(j = 1:ncol(dat_impute)), func(j))
save.image("../results/Writeup_revision10_main_gene_continuum_analysis3.RData")


# ###
# xx <- order(sapply(1:ncol(dat_impute), function(j){
#   mean(dat_ordered2[idx_trajectory2,j]) - mean(dat_ordered1[idx_trajectory1,j])
#   # mean(dat_ordered1[idx_trajectory1,j]) - mean(dat_ordered2[idx_trajectory2,j])
# }), decreasing = T)
# xx[1:5]
#
# par(mfrow = c(1,2))
# j <- 375
# col_vec <- rep("black", nrow(dat_ordered1)); col_vec[idx_trajectory1] <- "red"
# plot(x1, dat_ordered1[,j], col = col_vec, pch = 16, cex = 0.5)
# col_vec <- rep("black", nrow(dat_ordered2)); col_vec[idx_trajectory2] <- "red"
# plot(x2, dat_ordered2[,j], col = col_vec, pch = 16, cex = 0.5)
#
#
