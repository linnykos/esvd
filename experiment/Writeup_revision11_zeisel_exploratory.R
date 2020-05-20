rm(list=ls())
load("../results/step1_zeisel_gaussian_fitting.RData")

rgl::plot3d(svd_embedding[,1], svd_embedding[,2], svd_embedding[,3], asp = T,
            col = as.numeric(label_vec))

mat <- svd_embedding
cluster_labels <- as.numeric(label_vec)
neigh_size <- determine_minimium_neighborhood_size(mat)
neigh_size
set.seed(10)
zz <- compute_purity(mat, cluster_labels, neigh_size, num_samples = 500)
zz$avg_val
