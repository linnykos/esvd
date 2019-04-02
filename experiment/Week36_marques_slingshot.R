rm(list=ls())
load("../results/step4_factorization_spca.RData")
zz <- svd(dat)
d <- 3
naive <- zz$u[,1:d]%*%diag(zz$d[1:d])
u_mat <- res_our$u_mat[,1:d]

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

upscale_vec <- rep(NA, length(unique(cluster_labels)))
size_vec <- sapply(cluster_group_list, function(x){length(which(cluster_labels %in% x))})
for(i in 1:length(cluster_group_list)){
  upscale_vec[cluster_group_list[[i]]] <- (max(size_vec)/size_vec[i])^(1/2)
}

slingshot_func <- function(dat){
  reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*.25
  singlecell::slingshot(dat/reduction_factor, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
            cluster_group_list = cluster_group_list, verbose = T, b = 1,
            upscale_vec = upscale_vec)
}

naive_curves <- slingshot_func(naive)
our_curves <- slingshot_func(u_mat)

save.image("Week36_marques_slingshot.RData")


