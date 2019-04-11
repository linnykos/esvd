rm(list=ls())
load("../results/step5_clustering_spca.RData")
#
# naive_tube_list <- lapply(1:length(naive_curves$curves), function(x){
#   s_mat <- naive_curves$curves[[x]]$s[naive_curves$curves[[x]]$ord,]
#   construct_3d_tube(s_mat, radius = naive_curves$sd_val)
# })

x = 1
s_mat <- naive_curves$curves[[x]]$s[naive_curves$curves[[x]]$ord,]
radius = naive_curves$sd_val
