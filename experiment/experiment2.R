rm(list=ls())

set.seed(10)
library(eSVD)

suffix <- ""
ncores <- 20
doMC::registerDoMC(cores = ncores)

load(paste0("../results/step5_trajectory", suffix, ".RData"))

####

dat <- esvd_embedding$u_mat[,1:p]
trials = 100
reduction_percentage = 0.25
cores = ncores
verbose = T
starting_cluster = cluster_group_list[[1]][1]

x = 70
set.seed(10*x)
if(verbose & x %% floor(trials/10) == 0) print('*')
dat2 <- dat
for(i in 1:length(unique(cluster_labels))){
  idx <- which(cluster_labels == i)
  idx2 <- sample(idx, length(idx), replace = T)
  dat2[idx,] <- dat[idx2,]
}

zz <- eSVD::slingshot(dat2, cluster_labels, starting_cluster = starting_cluster,
          cluster_group_list = cluster_group_list, reduction_percentage = reduction_percentage,
          verbose = T)

#########

dat2 <- dat
reduction_percentage = 0.1
use_initialization = F
shrink = 1
thresh = 0.001
max_iter = 15
upscale_vec = NA

lineages <- eSVD:::.get_lineages(dat, cluster_labels, starting_cluster = starting_cluster,
                          cluster_group_list = cluster_group_list,
                          use_initialization = use_initialization)

reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*reduction_percentage
dat2 <- dat/reduction_factor

res <- eSVD:::.get_curves(dat2, cluster_labels, lineages, shrink = shrink,
                   thresh = thresh, max_iter = max_iter, upscale_vec = upscale_vec,
                   verbose = verbose)

curves <- res$pcurve_list

# adjust up
for(k in 1:length(curves)){
  curves[[k]]$s <- curves[[k]]$s*reduction_factor
}
