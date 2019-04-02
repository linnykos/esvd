rm(list=ls())
load("../results/step4_factorization_spca.RData")

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

d <- 3
dat <- res_our$u_mat[,1:d]
reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*.2
set.seed(10)
# curves <- slingshot(dat/reduction_factor, cluster_labels, starting_cluster = cluster_group_list[[1]][1], cluster_group_list = cluster_group_list, verbose = T,
#                     b = 1)

############

dat = dat/reduction_factor
shrink = 1
thresh = 0.001
max_iter = 15
b = 1
upscale_vec = NA
verbose = F
