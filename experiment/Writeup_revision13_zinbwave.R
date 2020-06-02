rm(list=ls())
load("../results/step7_additional_analyses_original.RData")

set.seed(10)
zinbwave_curves <- eSVD::slingshot(zinbwave_embedding[,1:3], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                                     cluster_group_list = cluster_group_list, shrink = 2,
                                     verbose = T, upscale_factor = 1, stretch = 2, max_iter = 3,
                                     squared = F)
