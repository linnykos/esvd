load(paste0("../results/step4_factorization", suffix, ".RData"))

cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

upscale_factor <- 1
reduction_percentage <- 0.2

p <- 3
set.seed(10)
esvd_curves <- eSVD::slingshot(esvd_embedding$u_mat[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                               cluster_group_list = cluster_group_list,
                               verbose = T, upscale_factor = upscale_factor,
                               reduction_percentage = reduction_percentage,
                               squared = T)

print(paste0(Sys.time(), ": Finished eSVD trajectory"))
save.image(paste0("../results/step5_trajectory", suffix, ".RData"))

# set.seed(10)
# esvd_bootstrap_list <- eSVD::bootstrap_curves(esvd_embedding$u_mat[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
#                                               cluster_group_list = cluster_group_list, trials = 100,
#                                               upscale_factor = upscale_factor, reduction_percentage = reduction_percentage, cores = ncores,
#                                               verbose = T, squared = F)
#
# print(paste0(Sys.time(), ": Finished eSVD bootstrap"))
# save.image(paste0("../results/step5_trajectory", suffix, ".RData"))
#
# esvd_sd_val <- eSVD::compute_curve_sd(esvd_curves, esvd_bootstrap_list, cores = ncores, verbose = T)
#
# print(paste0(Sys.time(), ": Finished eSVD standard val"))
# save.image(paste0("../results/step5_trajectory", suffix, ".RData"))

#########

set.seed(10)
svd_curves <- slingshot(svd_embedding[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                        cluster_group_list = cluster_group_list,
                        verbose = T, upscale_factor = upscale_factor,
                        reduction_percentage = reduction_percentage,
                        squared = T)

print(paste0(Sys.time(), ": Finished SVD trajectory"))
save.image(paste0("../results/step5_trajectory", suffix, ".RData"))

# set.seed(10)
# svd_bootstrap_list <- eSVD::bootstrap_curves(svd_embedding[,1:p], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
#                                              cluster_group_list = cluster_group_list, trials = 100,
#                                              upscale_factor = upscale_factor, reduction_percentage = reduction_percentage, cores = ncores,
#                                              verbose = T, squared = F)
#
#
# print(paste0(Sys.time(), ": Finished SVD bootstrap"))
# save.image(paste0("../results/step5_trajectory", suffix, ".RData"))
#
# svd_sd_val <- eSVD::compute_curve_sd(svd_curves, svd_bootstrap_list, cores = ncores, verbose = T)
#
# print(paste0(Sys.time(), ": Finished SVD standard val"))
# save.image(paste0("../results/step5_trajectory", suffix, ".RData"))

rm(list = c("tmp"))
print(paste0(Sys.time(), ": Finished trajectory"))
source_code_info <- c(source_code_info, readLines("../main/step5_trajectory.R"))
save.image(paste0("../results/step5_trajectory", suffix, ".RData"))
print(warnings())

