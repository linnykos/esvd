load(paste0("../results/step4_factorization", suffix, ".RData"))
session_info <- sessionInfo(); date_of_run <- Sys.time()

cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

set.seed(10)
esvd_curves_long <- eSVD::slingshot(esvd_embedding$u_mat, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                               cluster_group_list = cluster_group_list, shrink = 3,
                               verbose = T, upscale_factor = 1, stretch = 9999, max_iter = 3,
                               squared = T)

set.seed(10)
esvd_curves_short <- eSVD::slingshot(esvd_embedding$u_mat, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                                    cluster_group_list = cluster_group_list, shrink = 2,
                                    verbose = T, upscale_factor = 1, stretch = 2, max_iter = 3,
                                    squared = T)

print(paste0(Sys.time(), ": Finished eSVD trajectory"))
save.image(paste0("../results/step5_trajectory", suffix, ".RData"))

set.seed(10)
esvd_bootstrap_list <- eSVD::bootstrap_curves(esvd_embedding$u_mat, cluster_labels, lineages = esvd_curves_short$lineages,
                                              trials = 100, shrink = 2, upscale_factor = 1, stretch = 2, max_iter = 3,
                                              ncores = ncores, verbose = T)

print(paste0(Sys.time(), ": Finished eSVD bootstrap"))
save.image(paste0("../results/step5_trajectory", suffix, ".RData"))

esvd_sd_res <- eSVD::compute_curve_sd(esvd_curves_short$curves, esvd_bootstrap_list,
                                      ncores = ncores, verbose = T)
esvd_width <- max(sapply(1:2, function(i){
  quantile(apply(esvd_sd_res$mat_list[[i]], 1, function(x){quantile(x, probs = 0.95)}), probs = 0.95)
}))

print(paste0(Sys.time(), ": Finished eSVD standard val"))
save.image(paste0("../results/step5_trajectory", suffix, ".RData"))

#########

set.seed(10)
svd_curves_short <- slingshot(svd_embedding, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                        cluster_group_list = cluster_group_list,
                        verbose = T, upscale_factor = 1, shrink = 3, stretch = 2, max_iter = 3,
                        squared = T)

print(paste0(Sys.time(), ": Finished SVD trajectory"))
save.image(paste0("../results/step5_trajectory", suffix, ".RData"))

set.seed(10)
svd_bootstrap_list <- eSVD::bootstrap_curves(svd_embedding, cluster_labels, lineages = svd_curves_short$lineages,
                                             trials = 100, upscale_factor = 1, shrink = 3, stretch = 2, max_iter = 3,
                                             ncores = ncores, verbose = T)


print(paste0(Sys.time(), ": Finished SVD bootstrap"))
save.image(paste0("../results/step5_trajectory", suffix, ".RData"))

svd_sd_res <- eSVD::compute_curve_sd(svd_curves_short$curves, svd_bootstrap_list,
                                     ncores = ncores, verbose = T)

svd_width <- max(sapply(1:2, function(i){
  quantile(apply(svd_sd_res$mat_list[[i]], 1, function(x){quantile(x, probs = 0.95)})[200:900], probs = 0.95)
}))

print(paste0(Sys.time(), ": Finished SVD standard val"))
save.image(paste0("../results/step5_trajectory", suffix, ".RData"))

rm(list = c("tmp"))
print(paste0(Sys.time(), ": Finished trajectory"))
source_code_info <- c(source_code_info, readLines("../main/step5_trajectory.R"))
save.image(paste0("../results/step5_trajectory", suffix, ".RData"))
print(warnings())

