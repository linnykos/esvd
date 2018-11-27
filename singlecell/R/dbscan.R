dbscan <- function(dat, neighbor_count = 10, upper_cutoff = 14,
                   size_cutoff = 19){
  dist_mat <- as.matrix(dist(dat))
  min_dist <- stats::quantile(apply(dist_mat, 1, function(x){
    stats::quantile(x, 0.05)
  }), 0.75)
  n <- nrow(dat)

  adj_mat <- matrix(0, n, n)
  for(i in 1:n){
    idx1 <- which(dist_mat[,i] <= min_dist)
    idx2 <- order(dist_mat[,i], decreasing = F)[1:neighbor_count]
    idx <- intersect(idx1, idx2)
    if(length(idx) > 0) adj_mat[idx,i] <- 1
  }
  adj_mat <- adj_mat + t(adj_mat)
  adj_mat[which(adj_mat > 0)] <- 1

  # create the graph
  neigh_vec <- apply(adj_mat, 2, function(x){length(which(x != 0))})
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
  g <- igraph::simplify(g)
  assign_vec <- rep(NA, n)

  count_idx <- 0
  for(i in 1:n){
    if(!is.na(assign_vec[i])) next()
    if(neigh_vec[i] < upper_cutoff){assign_vec[i] <- 0; next()}

    q <- dequer::queue()
    count_idx <- count_idx + 1
    assign_vec[i] <- count_idx
    tmp <- sapply(as.numeric(igraph::neighbors(g, i)), function(x){
      dequer::pushback(q, x); invisible()
    })

    while(length(q) > 0){
      item <- dequer::pop(q)
      if(!is.na(assign_vec[item]) & assign_vec[item] != 0) next()
      assign_vec[item] <- count_idx
      if(neigh_vec[item] >= upper_cutoff){
        tmp <- sapply(as.numeric(igraph::neighbors(g, item)), function(x){
          dequer::pushback(q, x); invisible()
        })
      }
    }
  }

  # remove small clusters
  tab <- table(assign_vec)
  cluster_remove <- as.numeric(names(tab)[which(tab <= size_cutoff)])
  assign_vec[assign_vec %in% cluster_remove] <- 0
  mat <- cbind(sort(unique(assign_vec))[-1], 1:(length(unique(assign_vec))-1))
  for(i in 1:nrow(mat)){
    assign_vec[which(assign_vec == mat[i,1])] <- mat[i,2]
  }

  assign_vec[assign_vec == 0] <- NA
  assign_vec
}
