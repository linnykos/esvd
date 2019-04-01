rm(list=ls())
load("../results/step4_factorization_spca.RData")
zz <- svd(dat)
naive <- zz$u[,1:3]%*%diag(zz$d[1:3])
dat <- naive

cell_type_vec <- as.character(marques$cell.info$cell.type[cell_idx])
cell_type_vec <- as.factor(cell_type_vec)
cluster_labels <- as.numeric(cell_type_vec)
order_vec <- c("PP", "OP", "CO", "NF", "MF", "MO")
cluster_group_list <- lapply(order_vec, function(x){
  grep(paste0("^", x), levels(cell_type_vec))
})

##########

func <- function(x){
  set.seed(10*x)
  dat2 <- dat
  for(i in 1:length(unique(cluster_labels))){
    idx <- which(cluster_labels == i)
    idx2 <- sample(idx, length(idx), replace = T)
    dat2[idx,] <- dat[idx2,]
  }

  singlecell::slingshot(dat, cluster_labels, starting_cluster = cluster_group_list[[1]][1],
            cluster_group_list = cluster_group_list, verbose = F,
            b = max(apply(dat, 2, function(x){diff(range(x))}))/10)
}

# run in parallel
trials <- 100
res_list_naive <- foreach::"%dopar%"(foreach::foreach(x = 1:trials), func(x))

save.image("tmp2.RData")
save.image("../experiment/Week36_slingshot_bootstrap.RData")

##############

dat <- res_our$u_mat[,1:3]
res_list_our <- foreach::"%dopar%"(foreach::foreach(x = 1:trials), func(x))

save.image("tmp2.RData")
save.image("../experiment/Week36_slingshot_bootstrap.RData")


