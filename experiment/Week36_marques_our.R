rm(list=ls())
load("../results/step4_factorization_spca.RData")

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

u_mat <- res_our$u_mat

lineages <- .get_lineages(u_mat, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                          cluster_group_list = cluster_group_list)

# let's try bootstrapping a bit
trials <- 100
curve_list <- lapply(1:trials, function(x){
  print(x)
  set.seed(10*x)
  u_mat2 <- u_mat
  for(i in 1:length(unique(cluster_labels))){
    idx <- which(cluster_labels == i)
    idx2 <- sample(idx, length(idx), replace = T)
    u_mat2[idx,] <- u_mat[idx2,]
  }

  .get_lineages(u_mat2, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                cluster_group_list = cluster_group_list)

})

# we need a way to get an confidence tube hm...
