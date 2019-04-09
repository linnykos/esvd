rm(list=ls())
load("../results/step4_factorization_spca.RData")

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

our_curves <- singlecell::slingshot(res_our$u_mat[,1:3], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                                cluster_group_list = cluster_group_list,
                                verbose = T)

our_bootstrap_list <- singlecell::bootstrap_curves(res_our$u_mat[,1:3], cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                                   cluster_group_list = cluster_group_list, trials = 100)

our_sd_val <- singlecell::compute_curve_sd(our_curves, our_bootstrap_list)
