.initial_radius <- function(dat){
  dist_mat <- as.matrix(stats::dist(dat))
  max(apply(dist_mat, 2, function(x){min(x[x > 0])}))
}

.construct_graph <- function(dat, radius){
  n <- nrow(dat)
  dist_mat <- as.matrix(stats::dist(dat))
  adj_mat <- matrix(0, n, n)

  for(i in 1:n){
    idx <- which(dist_mat[i,] <= radius)
    if(length(idx) > 0){
      adj_mat[i,idx] <- 1; adj_mat[idx,i] <- 1
    }
  }

  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
  igraph::simplify(g)
}

# add the minimum amount of edges to get a connected graph
.adjust_graph <- function(dat, g){
  g_new <- g
  dist_mat <- as.matrix(stats::dist(dat))

  while(TRUE){
    comp_res <- igraph::components(g_new)
    k <- max(comp_res$membership)
    if(k == 1) break()

    # find which two clusters to merge
    combn_mat <- utils::combn(k, 2)
    dist_vec <- sapply(1:ncol(combn_mat), function(x){
      min(dist_mat[which(comp_res$membership == combn_mat[1,x]),
                   which(comp_res$membership == combn_mat[2,x])])
    })

    # find the specific pair of points to merge
    idx <- which.min(dist_vec)
    set1 <- which(comp_res$membership == combn_mat[1,idx])
    set2 <- which(comp_res$membership == combn_mat[2,idx])
    idx_pair <- which(dist_mat[set1, set2] == min(dist_vec), arr.ind = T)
    point_pair <- c(set1[idx_pair[1]], set2[idx_pair[2]])

    # do the merge
    g_new <- igraph::add_edges(g_new, point_pair)
  }

  g_new
}

.purity_point_pair <- function(g, label_vec, point_1, point_2){
  stopifnot(label_vec[point_1] == label_vec[point_2])

  path <- as.numeric(igraph::shortest_paths(g, from = point_1, to = point_2)$vpath[[1]])
  length(which(label_vec[path] == label_vec[point_1]))/length(label_vec)
}

.purity_data <- function(g, label_vec){
  k <- length(unique(label_vec))

  purity_list <- lapply(1:k, function(x){
    idx <- which(label_vec == x)
    combn_mat <- utils::combn(length(idx), 2)

    purity_vec <- sapply(1:ncol(combn_mat), function(i){
      .purity_point_pair(g, label_vec, idx[combn_mat[1,i]], idx[combn_mat[2,i]])
    })
  })

  mean(sapply(purity_list, mean))
}
