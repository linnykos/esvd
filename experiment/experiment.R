rm(list=ls())
load("../results/step5_trajectory_original.RData")

upscale_factor <- 1
reduction_percentage <- 0.2

p <- 3
set.seed(10)
esvd_curves <- eSVD::slingshot(esvd_embedding$u_mat[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                               cluster_group_list = cluster_group_list,
                               verbose = T, upscale_factor = upscale_factor,
                               reduction_percentage = reduction_percentage,
                               squared = T)

# checking of all points are accounted for
idx <- which(cluster_labels %in% esvd_curves$lineages[[1]])
## first extract which indices are in the first lineage
idx2 <- which(esvd_curves$curves[[1]]$W != 0)
all(idx %in% esvd_curves$idx[idx2])

## same thing for second lineage
idx <- which(cluster_labels %in% esvd_curves$lineages[[2]])
idx2 <- which(esvd_curves$curves[[2]]$W != 0)
all(idx %in% esvd_curves$idx[idx2])

###################################

slingshot_res <- esvd_curves
pseudotime_df_list <- .extract_pseudotimes(slingshot_res)
shared_df <- .compile_common_cells(pseudotime_df_list)

table(shared_df$consensus, cluster_labels[shared_df$cell_idx])

