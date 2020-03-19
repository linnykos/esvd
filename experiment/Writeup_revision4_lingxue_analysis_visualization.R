rm(list=ls())
load("../results/lingxue_analysis.RData")

# interactive 3D plot
# https://rpubs.com/aagarwal29/179912

i <- 7 # which dataset?
j <- 1 # how many dimensions?
k <- 2 # which distribution?

u_mat <- fit_all_list[[i]][[j]][[k]]$u_mat

rgl::plot3d(u_mat[,1], u_mat[,2], u_mat[,3], asp = T, col = as.numeric(preprocessing_list[[i]]$label_vec))
