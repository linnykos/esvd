load(paste0("../results/step5_trajectory", suffix, ".RData"))
session_info <- sessionInfo(); date_of_run <- Sys.time()

min_traj_pseudotime <- 6
segmentation_prep <- eSVD::prepare_data_for_segmentation(dat_impute, cluster_labels = cluster_labels,
                                                         curve_list = esvd_curves_short, min_traj_pseudotime = min_traj_pseudotime,
                                                         cluster_removal_idx_vec = c(2,3,4), cluster_removal_time_vec = rep(4,3))

segmentation_res <- eSVD::segment_genes_along_trajectories(segmentation_prep$dat1, segmentation_prep$dat2,
                                                           common_n = length(cell_idx_common$cell_idx_common),
                                                           standardize = T, verbose = T, ncores = 20)

print(paste0(Sys.time(), ": Finished cascading results"))
source_code_info <- c(source_code_info, readLines("../main/step6_cascade.R"))
save.image(paste0("../results/step6_cascade", suffix, ".RData"))
print(warnings())

