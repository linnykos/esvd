rm(list=ls())
load("../results/step3_zeisel_factorization.RData")

# u_mat <- svd_embedding
# u_mat <- esvd_missing_list[[11]][[1]]$u_mat
u_mat <- esvd_embedding$u_mat
rgl::plot3d(u_mat[,1], u_mat[,2], u_mat[,3], asp = T, col = as.numeric(label_vec),
            main = "Zeisel")

######3

mat <- esvd_missing_list[[15]][[1]]$u_mat
cluster_labels <- as.numeric(label_vec)
neigh_size <- determine_minimium_neighborhood_size(mat)
print(neigh_size)
set.seed(10)
zz <- compute_purity(mat, cluster_labels, neigh_size, num_samples = 1000)
print(zz$avg_val)

mat <- svd_embedding
cluster_labels <- as.numeric(label_vec)
neigh_size <- determine_minimium_neighborhood_size(mat)
print(neigh_size)
set.seed(10)
zz <- compute_purity(mat, cluster_labels, neigh_size, num_samples = 1000)
print(zz$avg_val)


########
idx <- c(7, 5, 12, 16, 20, 21)

for(i in idx){
  print(i)
  set.seed(10)

  mat <- esvd_missing_list[[i]][[1]]$u_mat
  cluster_labels <- as.numeric(label_vec)
  neigh_size <- determine_minimium_neighborhood_size(mat)
  print(neigh_size)
  set.seed(10)
  zz <- compute_purity(mat, cluster_labels, neigh_size)
  print(zz$avg_val)
}

for(i in c(1:3)){
  print(i)
  set.seed(10)

  mat <- esvd_missing_list[[12]][[i]]$u_mat
  cluster_labels <- as.numeric(label_vec)
  neigh_size <- determine_minimium_neighborhood_size(mat)
  print(neigh_size)
  set.seed(10)
  zz <- compute_purity(mat, cluster_labels, neigh_size)
  print(zz$avg_val)
}



mat <- svd_embedding
cluster_labels <- as.numeric(label_vec)
neigh_size <- determine_minimium_neighborhood_size(mat)
print(neigh_size)
set.seed(10)
zz <- compute_purity(mat, cluster_labels, neigh_size)
print(zz$avg_val)
