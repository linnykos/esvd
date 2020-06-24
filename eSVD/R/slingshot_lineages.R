#' Estimate the lineages (via Slingshot)
#'
#' Code adapted from https://github.com/kstreet13/slingshot. See
#' our paper for a list of changes to adapt the original Slingshot code
#' to our setting
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}. Cannot include \code{NA}
#' @param starting_cluster the "origin" cluster that all the lineages will start
#' from
#' @param cluster_group_list list denoting the hierarchy and order of the clusters
#' @param squared boolean on whether or not to square the distance matrix
#'
#' @return A list of cluster indices, with \code{starting_cluster} starting as
#' its first element
.get_lineages <- function(dat, cluster_labels, starting_cluster,
                          cluster_group_list = NA, squared = F){
  stopifnot(!is.list(cluster_group_list) || starting_cluster %in% cluster_group_list[[1]])
  stopifnot(all(cluster_labels > 0), all(cluster_labels %% 1 == 0), length(unique(cluster_labels)) == max(cluster_labels))
  if(all(!is.na(cluster_group_list))){
    tmp <- unlist(cluster_group_list)
    stopifnot(length(tmp) == length(unique(tmp)), length(tmp) == length(unique(cluster_labels)))
  }

  ### construct the distance matrix
  dist_mat <- .compute_cluster_distances(dat, cluster_labels)
  if(squared) dist_mat <- dist_mat^2

  if(all(is.na(cluster_group_list))){
    ### construct the spt (shortest path tree)
    g <- .construct_graph(dist_mat)

    ### identify lineages (paths through trees)
    lineages <- .construct_lineages(g, starting_cluster = starting_cluster)

  } else {
    lineages <- .construct_lineage_from_hierarchy(dist_mat, cluster_group_list,
                                                  starting_cluster)
  }

  lineages
}

#############

#' Compute the Hotelling t test statistic for non-equal variances
#'
#' See \url{http://www.real-statistics.com/multivariate-statistics/hotellings-t-square-statistic/hotellings-t-square-unequal-covariance-matrices/}
#'
#' @param mean_vec1 vector
#' @param cov_mat1 matrix
#' @param n1 numeric
#' @param mean_vec2 vector
#' @param cov_mat2 matrix
#' @param n2 numeric
#' @param tol numeric
#'
#' @return numeric
.covariance_distance <- function(mean_vec1, cov_mat1, n1, mean_vec2, cov_mat2, n2, tol = 1e-5){
  mat <- cov_mat1/n1 + cov_mat2/n2

  if(Matrix::rankMatrix(mat) < ncol(mat)){
    eigen_res <- eigen(mat)
    eigen_res$values[eigen_res$values < tol] <- tol
    inv_mat <- eigen_res$vectors %*% diag(1/eigen_res$values) %*% t(eigen_res$vectors)
  } else {
    inv_mat <- solve(mat)
  }

  as.numeric(t(mean_vec1 - mean_vec2) %*% inv_mat %*% (mean_vec1 - mean_vec2))
}

#' Compute distances between all the clusters
#'
#' @param dat matrix
#' @param cluster_labels vector with length equal to \code{nrow(dat)}
#'
#' @return matrix
.compute_cluster_distances <- function(dat, cluster_labels){
  k <- max(cluster_labels)
  dist_mat <- matrix(0, k, k)

  for(i in 1:(k-1)){
    for(j in (i+1):k){
      idx1 <- which(cluster_labels == i)
      idx2 <- which(cluster_labels == j)

      mean_vec1 <- colMeans(dat[idx1,])
      mean_vec2 <- colMeans(dat[idx2,])

      cov_mat1 <- stats::cov(dat[idx1,])
      cov_mat2 <- stats::cov(dat[idx2,])

      n1 <- length(idx1); n2 <- length(idx2)

      dist_mat[i,j] <- .covariance_distance(mean_vec1, cov_mat1, n1, mean_vec2, cov_mat2, n2)
      dist_mat[j,i] <- dist_mat[i,j]
    }
  }

  dist_mat
}

#' Construct lineages from hiearchy
#'
#' Using the order of clusters in \code{cluster_group_list}, construct the
#' lineage based on \code{dist_mat}. The function operates by appending clusters onto
#' existing trees. This operates in rounds equal to \code{length(cluster_group_list)}.
#' Specifically, within each round, the function first uses \code{eSVD:::.enumerate_dist_from_trees}
#' to populate a distance matrix from the previous round's trees, and then uses
#' \code{eSVD:::.enumerate_dist_between_levels} and \code{eSVD:::.enumerate_dist_within_levels} to
#' populate the new edges connecting the previous round's clusters to this round's clusters, as well
#' as this round's clusters to each other. Then, it constructs a list of trees via
#' \code{igraph::shortest_paths} and uses \code{eSVD:::.find_all_unique_paths} to "prune" the list
#' of trees by removing paths that are strictly contained in other paths. Then it proceeds to the next round.
#'
#' @param dist_mat (symmetric) distance matrix
#' @param cluster_group_list list
#' @param starting_cluster numeric
#'
#' @return list, enumerating the different lineages
.construct_lineage_from_hierarchy <- function(dist_mat, cluster_group_list,
               starting_cluster = 1){
  n <- max(unlist(cluster_group_list))
  len <- length(cluster_group_list)
  stopifnot(len > 1)
  tree_list <- list(c(starting_cluster))

  # iterate through all the levels in the cluster_group_list
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
      g <- igraph::add_edges(g, edges = c(edge_mat[j,1], edge_mat[j,2]),
                             attr = list(weight = edge_mat[j,3]))
    }

    # find shortest path tree
    path_list <- suppressWarnings(igraph::shortest_paths(g, from = starting_cluster,
                                        output = "vpath")$vpath)

    # find all unique paths
    tree_list <- .find_all_unique_paths(path_list, starting_cluster)
  }

  names(tree_list) <- paste('Lineage', seq_along(tree_list), sep='')

  tree_list
}

#' Enumerate distances from paths stored in a list
#'
#' @param dist_mat (symmetric) distance matrix
#' @param tree_list list of paths, corresponding to indices from 1 to \code{nrow(dist_mat)}
#'
#' @return a matrix with 3 columns, representing the two indicies and the value from \code{dist_mat}
.enumerate_dist_from_trees <- function(dist_mat, tree_list){
  if(all(sapply(tree_list, length) == 1)) return(numeric(0))

  edge_list <- lapply(1:length(tree_list), function(i){
    if(length(tree_list[[i]]) > 1){

      t(sapply(1:(length(tree_list[[i]])-1), function(j){
        idx1 <- tree_list[[i]][j]
        idx2 <- tree_list[[i]][j+1]

        c(idx1, idx2, dist_mat[idx1, idx2])
      }))
    }
  })

  do.call(rbind, edge_list)
}

#' Enumerate distances between paths and a set of indices
#'
#' Enumerate edges corresponding to the last index in each element in \code{tree_list}
#' and any index in \code{cluster_vec}
#'
#' @param dist_mat (symmetric) distance matrix
#' @param tree_list list of paths, corresponding to indices from 1 to \code{nrow(dist_mat)}
#' @param cluster_vec list of indices from 1 to \code{nrow(dist_mat)}
#'
#' @return a matrix with 3 columns, representing the two indicies and the value from \code{dist_mat}
.enumerate_dist_between_levels <- function(dist_mat, tree_list, cluster_vec){
  # enumerate all the leaves, one for each tree
  leaf_vec <- sapply(tree_list, function(x){x[length(x)]})
  stopifnot(length(leaf_vec) == length(unique(leaf_vec)))

  node_pairings <- as.matrix(expand.grid(leaf_vec, cluster_vec))

  edge_mat <- t(sapply(1:nrow(node_pairings), function(i){
    idx1 <- node_pairings[i,1]
    idx2 <- node_pairings[i,2]

    c(idx1, idx2, dist_mat[idx1, idx2])
  }))

  edge_mat
}

#' Enumerate distances within a vector of indices
#'
#' @param dist_mat (symmetric) distance matrix
#' @param cluster_vec list of indices from 1 to \code{nrow(dist_mat)}
#'
#' @return a matrix with 3 columns, representing the two indicies and the value from \code{dist_mat}
.enumerate_dist_within_levels <- function(dist_mat, cluster_vec){
  if(length(cluster_vec) < 2) return(numeric(0))
  combn_mat <- utils::combn(length(cluster_vec), 2)

  edge_mat <- lapply(1:ncol(combn_mat), function(i){
    idx1 <- cluster_vec[combn_mat[1,i]]
    idx2 <- cluster_vec[combn_mat[2,i]]

    matrix(c(idx1, idx2, dist_mat[idx1, idx2], idx2, idx1, dist_mat[idx1, idx2]),
           nrow = 2, ncol = 3, byrow = T)
  })

  do.call(rbind, edge_mat)
}

#' Create an undirected graph from a distance matrix
#'
#' @param dist_mat (symmetric) distance matrix
#'
#' @return \code{igraph} object
.construct_graph <- function(dist_mat){
  n <- nrow(dist_mat)
  combn_mat <- utils::combn(n, 2)
  g <- igraph::graph.empty(n = n, directed = F)

  for(j in 1:ncol(combn_mat)){
    idx1 <- combn_mat[1,j]; idx2 <- combn_mat[2,j]

    g <- igraph::add_edges(g, edges = c(idx1, idx2),
                           attr = list(weight = dist_mat[idx1, idx2]))
  }

  g
}

#' Construct the lineages
#'
#' Enumerates the shortest paths from \code{starting_cluster} to
#' each of the leaves in \code{g}
#'
#' @param g \code{igraph} object
#' @param starting_cluster positive integer
#'
#' @return a list
.construct_lineages <- function(g, starting_cluster){
  path_list <- igraph::shortest_paths(g, from = starting_cluster,
                                     output = "vpath")$vpath

  lineages <- .find_all_unique_paths(path_list, starting_cluster)

  for(i in 1:length(lineages)){
    lineages[[i]] <- as.numeric(lineages[[i]])
  }

  names(lineages) <- paste('Lineage' ,seq_along(lineages), sep='')

  lineages
}

#' Remove duplicate paths
#'
#' We only keep paths in \code{path_list} that 1) start from \code{starting_cluster} and
#' 2) are not strictly contained in any other path in our outputed list
#'
#' @param path_list list of vectors of indices
#' @param starting_cluster index
#'
#' @return a pruned version of \code{path_list}
.find_all_unique_paths <- function(path_list, starting_cluster){
  # remove all paths that do not have the correct starting_cluster
  path_list <- path_list[which(sapply(path_list, function(x){x[1] == starting_cluster}))]
  stopifnot(length(path_list) >= 1)

  len <- length(path_list)

  # remove any paths strictly contained in another
  bool_vec <- sapply(1:len, function(x){
    !any(sapply(path_list[-x], function(y){
      all(path_list[[x]] %in% y)
    }))
  })

  path_list <- path_list[which(bool_vec)]

  lapply(path_list, function(x){
    as.numeric(x)
  })
}
