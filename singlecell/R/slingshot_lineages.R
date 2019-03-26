#' Estimate the lineages (via Slingshot)
#'
#' Note: I removed the functionality to explicitly label a ending cluster
#' (might put back in later?).
#'
#' Note: I removed the Omega parameter, which, to my understanding,
#' controls if a lineage (i.e. tree) is split into two separate lineages.
#' Currently, the only way for a forest to occur is if the KNN graph is
#' naturally disconnected.
#'
#' Note: Currently, this code is hard-coded to work for only a fully connected
#' KNN graph.
#'
#' Code adapted from https://github.com/kstreet13/slingshot.
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels, na.rm = T)}. Can include \code{NA}
#' @param starting_cluster the "origin" cluster that all the lineages will start
#' from
#' @param knn positive integer, possibly \code{NA}
#' @param remove_outlier boolean
#' @param percentage percentage dictating which points are considered outliers
#'
#' @return A list of cluster indices, with \code{starting_cluster} starting as
#' its first element
.get_lineages <- function(dat, cluster_labels, starting_cluster, cluster_group_list = NA){
  stopifnot(!is.list(cluster_group_list) || starting_cluster %in% cluster_group_list[[1]])

  ### construct the distance matrix
  dist_mat <- .compute_cluster_distances(dat, cluster_labels)

  ### construct the spt
  g <- .construct_graph_hierarchy(dist_mat, cluster_labels = cluster_labels,
                                  cluster_group_list = cluster_group_list)

  ### identify lineages (paths through trees)
  lineages <- .construct_lineages(g, starting_cluster = starting_cluster)

  lineages
}

#############

.covariance_distance <- function(mean_vec1, cov_mat1, mean_vec2, cov_mat2){
  as.numeric(t(mean_vec1 - mean_vec2) %*% solve(cov_mat1 + cov_mat2) %*% (mean_vec1 - mean_vec2))
}

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

      dist_mat[i,j] <- .covariance_distance(mean_vec1, cov_mat1, mean_vec2, cov_mat2)
      dist_mat[j,i] <- dist_mat[i,j]
    }
  }

  dist_mat
}

.construct_graph_hierarchy <- function(dist_mat, cluster_labels, cluster_group_list,
               starting_cluster = 1){

  n <- nrow(dist_mat)
  k <- length(cluster_group_list)
  g <- igraph::graph.empty(n = n, directed = T)

  # add edges within each group
  for(i in 1:k){
    m <- length(cluster_group_list[[i]])

    if(m > 1){
      combn_mat <- utils::combn(m, 2)
      for(j in 1:ncol(combn_mat)){
        idx1 <- cluster_group_list[[i]][combn_mat[1,j]]
        idx2 <- cluster_group_list[[i]][combn_mat[2,j]]

        g <- igraph::add_edges(g, edges = matrix(c(idx1, idx2, idx2, idx1), 2, 2),
                               attr = list(weight = dist_mat[idx1, idx2]))
      }
    }
  }

  # add edges between groups
  for(i in 1:(k-1)){
    edge_mat <- as.matrix(expand.grid(cluster_group_list[[i]], cluster_group_list[[i+1]]))
    for(j in 1:nrow(edge_mat)){
      idx1 <- edge_mat[j,1]
      idx2 <- edge_mat[j,2]

      g <- igraph::add_edges(g, edges = c(idx1, idx2),
                             attr = list(weight = dist_mat[idx1, idx2]))
    }
  }

  #cbind(igraph::as_edgelist(g, names = T), igraph::edge_attr(g, "weight", index = E(g)))
  g
}

#' Construct the lineages
#'
#' Enumerates the shortest paths from \code{starting_cluster} to
#' each of the leaves in \code{spt_graph}
#'
#' @param spt_graph \code{igraph} object
#' @param starting_cluster positive integer
#'
#' @return a list
.construct_lineages <- function(g, starting_cluster){
  path_list <- igraph::shortest_paths(g, from = starting_cluster,
                                     output = "vpath")$vpath

  lineages <- .find_all_unique_paths(path_list)

  for(i in 1:length(lineages)){
    lineages[[i]] <- as.numeric(lineages[[i]])
  }

  names(lineages) <- paste('Lineage' ,seq_along(lineages), sep='')

  lineages
}

.find_all_unique_paths <- function(path_list){
  len <- length(path_list)

  bool_vec <- sapply(1:len, function(x){
    !any(sapply(path_list[-x], function(y){
      all(path_list[[x]] %in% y)
    }))
  })

  path_list[which(bool_vec)]
}

.compute_ellipse_points <- function(mean_vec, cov_mat, scale_factor = 1){
  eig <- eigen(cov_mat)
  alpha <- atan(eig$vectors[2,1]/eig$vectors[1,1])
  if(alpha < 0) alpha <- alpha + 2*pi

  a <- sqrt(eig$values[1])*scale_factor
  b <- sqrt(eig$values[2])*scale_factor

  theta_grid <- seq(0, 2*pi, length.out = 100)
  ellipse_x <- a*cos(theta_grid)
  ellipse_y <- b*sin(theta_grid)

  R <- matrix(c(cos(alpha), -sin(alpha), sin(alpha), cos(alpha)), 2, 2)
  val <- cbind(ellipse_x, ellipse_y) %*% R

  val <- t(apply(val, 1, function(x){x + mean_vec}))
}


