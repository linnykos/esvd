load(paste0("../results/step4_factorization", suffix, ".RData"))

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

curves <- singlecell::slingshot(res_our$u_mat, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                                cluster_group_list = cluster_group_list)

rm(list = c("cell_type_vec"))
print(paste0(Sys.time(), ": Finished clustering"))
save.image(paste0("../results/step5_clustering", suffix, ".RData"))

