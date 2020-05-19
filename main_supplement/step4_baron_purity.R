set.seed(10)
load(paste0("../results/step3_baron_factorization", suffix, ".RData"))
session_info <- sessionInfo(); date_of_run <- Sys.time()

# compute the desired neighborhood sizes
neighborhood_size_vec <- sapply(1:length(preprocessing_list), function(i){
  print(i)
  mat <- svd_embedding_list[[i]]
  val1 <- eSVD:::determine_minimium_neighborhood_size(mat, verbose = F)

  mat <- esvd_embedding_list[[i]]$u_mat
  val2 <- eSVD:::determine_minimium_neighborhood_size(mat, verbose = F)

  max(val1, val2)
})

# compute the purity for each of the SVD embeddings
svd_purity <- sapply(1:length(preprocessing_list), function(i){
  set.seed(10)
  print(i)
  mat <- svd_embedding_list[[i]]
  cluster_labels <- as.numeric(preprocessing_list[[i]]$label_vec)
  purity <- eSVD:::compute_purity(mat, cluster_labels, neighborhood_size = neighborhood_size_vec[[i]])
  purity$avg_val
})

# compute the purity for each of the eSVD embeddings
esvd_purity <- sapply(1:length(preprocessing_list), function(i){
  set.seed(10)
  print(i)
  mat <- esvd_embedding_list[[i]]$u_mat
  cluster_labels <- as.numeric(preprocessing_list[[i]]$label_vec)
  purity <- eSVD:::compute_purity(mat, cluster_labels, neighborhood_size = neighborhood_size_vec[[i]])
  purity$avg_val
})
