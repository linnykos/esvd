compute_purity <- function(mat, cluster_labels, neighborhood_size, num_samples = 200){
  g <- .construct_neighborhood_graph(mat, neighborhood_size = neighborhood_size)
  stopifnot(igraph::components(g)$no == 1)

  k <- max(cluster_labels)
  idx_list <- lapply(1:k, function(i){which(cluster_labels == i)})

  value_list <- lapply(1:k, function(i){
    ni <- length(idx_list[[i]])
    combn_mat <- utils::combn(ni, 2)
    if(ncol(combn_mat) > num_samples){
      combn_mat <- combn_mat[,sample(1:ncol(combn_mat), num_samples)]
    }

    for(j in 1:2){combn_mat[j,] <- idx_list[[i]][combn_mat[j,]]}
    sapply(1:ncol(combn_mat), function(j){
      .compute_pairwise_purity(g, idx1 = combn_mat[1,j], idx2 = combn_mat[2,j], cluster_labels = cluster_labels)
    })
  })

  # list(avg_val = mean(sapply(value_list, mean)), value_list = value_list)
  list(avg_val = mean(unlist(value_list)), value_list = value_list)
}

determine_minimium_neighborhood_size <- function(mat, max_iter = 20, verbose = T){
  n <- nrow(mat)

  # binary search
  low <- 1
  high <- n-1

  g <- .construct_neighborhood_graph(mat, neighborhood_size = low)
  connected_low <- igraph::components(g)$no == 1
  if(connected_low) return(low)
  g <- .construct_neighborhood_graph(mat, neighborhood_size = high)
  connected_high <- igraph::components(g)$no == 1
  stopifnot(!connected_low & connected_high)

  iter <- 1
  while(iter < max_iter){
    mid <- round((low+high)/2)
    if(verbose) print(paste0("Iteration: ", iter, "//Low: ", low, "//Mid: ", mid, "//High: ", high))
    if(mid %in% c(low, high)) break()

    g <- .construct_neighborhood_graph(mat, neighborhood_size = mid)
    connected_mid <- igraph::components(g)$no == 1

    if(connected_mid == FALSE){
      low <- mid
      connected_low <- connected_mid
    } else {
      high <- mid
      connected_high <- connected_mid
    }

    iter <- iter+1
  }

  high
}

#############

.construct_neighborhood_graph <- function(mat, neighborhood_size){
  n <- nrow(mat)
  dist_mat <- as.matrix(stats::dist(mat))
  diag(dist_mat) <- Inf
  neighbor_list <- lapply(1:n, function(i){
    order(dist_mat[i,], decreasing = F)[1:neighborhood_size]
  })

  edge_mat <- do.call(cbind, lapply(1:n, function(i){
    rbind(rep(i, neighborhood_size), neighbor_list[[i]])
  }))

  g <- igraph::graph.empty(n = n, directed = F)
  g <- igraph::add_edges(g, edges = edge_mat)
  g <- igraph::simplify(g)
  g
}

.compute_pairwise_purity <- function(g, idx1, idx2, cluster_labels){
  stopifnot(cluster_labels[idx1] == cluster_labels[idx2])
  k <- cluster_labels[idx1]
  path <- igraph::shortest_paths(g, from = idx1, to = idx2)

  len <- length(path$vpath)
  val_vec <- sapply(1:len, function(i){
    vec <- path$vpath[[i]]
    if(length(vec) == 2) return(1)
    length(which(cluster_labels[vec[-c(1,length(vec))]] == k))/(length(vec)-2)
  })
  max(val_vec)
}
