rm(list=ls())
source("../experiment/Week27_simulation_generator.R")

set.seed(10)
simulation <- .data_generator(total = 100)
dat <- simulation$cell_mat
cluster_labels <- rep(1:simulation$h, each = simulation$n_each)

lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
cluster_mat <- res$cluster_mat
curves <- .get_curves_tmp(dat, cluster_mat, lineages, reassign = F)

#######

col_vec <- c(rgb(205,40,54,maxColorValue=255), #red
             rgb(180,200,255,maxColorValue=255), #purple
             rgb(100,100,200,maxColorValue=255), #blue
             rgb(149,219,144,maxColorValue=255)) #green
plot(simulation$cell_mat[,1], simulation$cell_mat[,2],
     xlim = range(c(simulation$cell_mat[,1], 0)),
     ylim = range(c(simulation$cell_mat[,2], 0)),
     col = col_vec[rep(1:simulation$h, each = simulation$n_each)], asp = T,
     pch = 16, xlab = "Latent dim. 1", ylab = "Latent dim. 2", main = "Cell latent vectors")


for(i in 1:length(curves)){
  ord <- curves[[i]]$ord
  lines(curves[[i]]$s[ord,1], curves[[i]]$s[ord,2])
}

# lines(s[,1], s[,2])

############################

shrink = 1; extend = 'y'; reweight = TRUE; reassign = F
thresh = 0.001; maxit = 15; stretch = 2
shrink_method = 'cosine'; allow.breaks = TRUE; b = 1

#####

cluster_mat <- .construct_cluster_matrix(cluster_labels)
k <- ncol(cluster_mat)
num_lineage <- length(grep("Lineage", names(lineages))) # number of lineages
clusters <- colnames(cluster_mat)
n <- nrow(dat)
p <- ncol(dat)
centers <- .compute_cluster_center(dat, cluster_mat)

W <- .initialize_weight_matrix(cluster_mat, lineages)
D <- matrix(NA, ncol = ncol(W), nrow = nrow(W))
