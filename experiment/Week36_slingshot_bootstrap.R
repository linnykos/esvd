rm(list=ls())
load("../results/step4_factorization_spca.RData")
zz <- svd(dat)
naive <- zz$u[,1:3]%*%diag(sqrt(zz$d[1:3]))
dat <- naive
reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*.25

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

##########

func <- function(x){
  reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*.25
  set.seed(10*x)
  dat2 <- dat
  for(i in 1:length(unique(cluster_labels))){
    idx <- which(cluster_labels == i)
    idx2 <- sample(idx, length(idx), replace = T)
    dat2[idx,] <- dat[idx2,]
  }

  singlecell::slingshot(dat2/reduction_factor, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
            cluster_group_list = cluster_group_list, verbose = F,
            b = 1, upscale_vec = upscale_vec)
}

# run in parallel
doMC::registerDoMC(cores = 15)
trials <- 50
res_list_naive <- foreach::"%dopar%"(foreach::foreach(x = 1:trials), func(x))

save.image("tmp2.RData")
save.image("../experiment/Week36_slingshot_bootstrap.RData")

##############

dat <- res_our$u_mat[,1:3]
res_list_our <- foreach::"%dopar%"(foreach::foreach(x = 1:trials), func(x))

save.image("tmp2.RData")
save.image("../experiment/Week36_slingshot_bootstrap.RData")


