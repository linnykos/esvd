rm(list=ls())
load("../results/step5_clustering.RData")

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

################

reduction_percentage <- 0.1
dat <- res_our$u_mat[,1:p]
reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*reduction_percentage
dat2 <- dat/reduction_factor

starting_cluster <- cluster_group_list[[1]][1]

###################

dat <- dat2
use_initialization = F
stopifnot(!is.list(cluster_group_list) || starting_cluster %in% cluster_group_list[[1]])
stopifnot(all(cluster_labels > 0), all(cluster_labels %% 1 == 0), length(unique(cluster_labels)) == max(cluster_labels))
if(all(!is.na(cluster_group_list))){
  tmp <- unlist(cluster_group_list)
  stopifnot(length(tmp) == length(unique(tmp)), length(tmp) == length(unique(cluster_labels)))
}

### construct the distance matrix
dist_mat <- .compute_cluster_distances(dat, cluster_labels)
if(use_initialization){
  bool_mat <- .initial_edges(dat, cluster_labels)
  if(all(rowSums(bool_mat) >= 1)){
    dist_mat[!bool_mat] <- Inf
  }
}

###########################

n <- length(unique(unlist(cluster_group_list)))
len <- length(cluster_group_list)
stopifnot(len > 1)
tree_list <- list(c(starting_cluster))

for(i in 2:len){
  # for each level k, first populate a new distance matrix with distances of all
  #  trees up until level k-1
  edge_mat1 <- .enumerate_dist_from_trees(dist_mat, tree_list)

  # populate distance matrix for all edges between one cluster in level k and
  #   one leaf cluster in level k-1
  stopifnot(!any(unlist(tree_list) %in% cluster_group_list[[i]]))
  edge_mat2 <- .enumerate_dist_between_levels(dist_mat, tree_list, cluster_group_list[[i]])

  # populate distance matrix between all clusters in level k
  edge_mat3 <- .enumerate_dist_within_levels(dist_mat, cluster_group_list[[i]])

  # populate igraph object
  edge_mat <- do.call(rbind, list(edge_mat1, edge_mat2, edge_mat3))
  g <- igraph::graph.empty(n = n, directed = T)

  for(j in 1:nrow(edge_mat)){
    idx1 <- edge_mat[j,1]; idx2 <- edge_mat[j,2]

    g <- igraph::add_edges(g, edges = c(idx1, idx2),
                           attr = list(weight = dist_mat[idx1, idx2]))
  }

  # find shortest path tree
  path_list <- suppressWarnings(igraph::shortest_paths(g, from = starting_cluster,
                                      output = "vpath")$vpath)

  # find all unique paths
  tree_list <- .find_all_unique_paths(path_list, starting_cluster)
}

tree_list
