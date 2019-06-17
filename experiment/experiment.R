rm(list=ls())
set.seed(10)

suffix <- ""
family <- "gaussian"
ncores <- 15
doMC::registerDoMC(cores = ncores)

load(paste0("../results/step4_factorization", suffix, ".RData"))

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

p <- 3
tmp <- svd(dat_impute)
naive_embedding <- tmp$u[,1:p] %*% diag(sqrt(tmp$d[1:p]))
naive_curves <- singlecell::slingshot(naive_embedding, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
                                      cluster_group_list = cluster_group_list,
                                      verbose = F, upscale_vec = upscale_vec)

##########

dat <- naive_curves

