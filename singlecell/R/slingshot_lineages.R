#' Estimate the lineages (via Slingshot)
#'
#' Note: I removed the functionality to explicitly label a ending cluster
#' (might put back in later?).
#'
#' Note: I removed the Omega parameter, which, to my understanding,
#' controls if a lineage (i.e. tree) is split into two separate lineages.
#'
#' Code adapted from https://github.com/kstreet13/slingshot.
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}. Cannot include \code{NA}
#' @param starting_cluster the "origin" cluster that all the lineages will start
#' from
#' @param cluster_group_list list denoting the hierarchy and order of the clusters
#'
#' @return A list of cluster indices, with \code{starting_cluster} starting as
#' its first element
.get_lineages <- function(dat, cluster_labels, starting_cluster, cluster_group_list = NA){
  stopifnot(!is.list(cluster_group_list) || starting_cluster %in% cluster_group_list[[1]])
  stopifnot(all(cluster_labels > 0), all(cluster_labels %% 1 == 0), length(unique(cluster_labels)) == max(cluster_labels))
  if(all(!is.na(cluster_group_list))){
    tmp <- unlist(cluster_group_list)
    stopifnot(length(tmp) == length(unique(tmp)), length(tmp) == length(unique(cluster_labels)))
  }

  ### construct the distance matrix
  dist_mat <- .compute_cluster_distances(dat, cluster_labels)

  ### construct the spt (shortest path tree)
  g <- .construct_graph_hierarchy(dist_mat, cluster_labels = cluster_labels,
                                  cluster_group_list = cluster_group_list)

  ### identify lineages (paths through trees)
  lineages <- .construct_lineages(g, starting_cluster = starting_cluster)

  lineages
}

#############

.covariance_distance <- function(mean_vec1, cov_mat1, mean_vec2, cov_mat2, tol = 0.1){
  mat <- cov_mat1 + cov_mat2

  eigen_res <- eigen(mat)
  eigen_res$values[eigen_res$values < tol] <- tol
  mat <- eigen_res$vectors %*% diag(1/eigen_res$values) %*% t(eigen_res$vectors)

  as.numeric(t(mean_vec1 - mean_vec2) %*% mat %*% (mean_vec1 - mean_vec2))
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

  g <- igraph::graph.empty(n = n, directed = T)

  edge_mat <- .populate_edge_matrix(cluster_labels, cluster_group_list)

  # populate the graph
  for(j in 1:ncol(edge_mat)){
    idx1 <- edge_mat[1,j]; idx2 <- edge_mat[2,j]

    g <- igraph::add_edges(g, edges = c(idx1, idx2),
                           attr = list(weight = dist_mat[idx1, idx2]))
  }

  g
}

.populate_edge_matrix <- function(cluster_labels, cluster_group_list = NA){
  if(all(!is.na(cluster_group_list))){
    k <- length(cluster_group_list)

    # add edges within each group
    edge_mat1 <- lapply(1:k, function(i){
      m <- length(cluster_group_list[[i]])
      if(m >= 2){
        combn_mat <- utils::combn(m, 2)
        combn_mat[1,] <- cluster_group_list[[i]][combn_mat[1,]]
        combn_mat[2,] <- cluster_group_list[[i]][combn_mat[2,]]
        combn_mat <- cbind(combn_mat, combn_mat[c(2,1),]) #the reverse edges
      } else {
        numeric(0)
      }
    })
    edge_mat1 <- do.call(cbind, edge_mat1)

    # add edges between groups
    edge_mat2 <- lapply(1:(k-1), function(i){
      t(as.matrix(expand.grid(cluster_group_list[[i]], cluster_group_list[[i+1]])))
    })
    edge_mat2 <- do.call(cbind, edge_mat2)
    edge_mat <- cbind(edge_mat1, edge_mat2)

  } else {
    edge_mat <- utils::combn(length(unique(cluster_labels)), 2)
    edge_mat <- cbind(edge_mat, edge_mat[c(2,1),]) #the reverse edges
  }

  edge_mat
}

#' Construct the lineages
#'
#' Enumerates the shortest paths from \code{starting_cluster} to
#' each of the leaves in \code{spt_graph}
#'
#' @param g \code{igraph} object
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

########

.binary_search <- function(dat, evaluation_func, target_val,
                           range_val = c(0,1), start_min = 0.1, start_max = 10,
                           max_iter = 10, ...){

  lower_val <- evaluation_func(dat, start_min, ...)
  upper_val <- evaluation_func(dat, start_max, ...)
  if(lower_val > target_val) return(start_min)
  if(upper_val < target_val) return(start_max)

  iter <- 1;
  min_val <- start_min; max_val <- start_max

  while(TRUE){
    try_val <- (min_val + max_val)/2
    mid_val <- evaluation_func(dat, try_val, ...)

    if(target_val < mid_val) {
      upper_val <- mid_val
    } else {
      lower_val <- mid_val
    }
  }

  try_val
}

.ellipse_points <- function(mean_vec, cov_mat){
  stopifnot(length(mean_vec) == 2, all(dim(cov_mat) == 2))

  eig <- eigen(cov_mat)
  alpha <- atan(eig$vectors[2,1]/eig$vectors[1,1])
  if(alpha < 0) alpha <- alpha + 2*pi

  a <- sqrt(eig$values[1])
  b <- sqrt(eig$values[2])

  theta_grid <- seq(0, 2*pi, length.out = 100)
  ellipse_x <- a*cos(theta_grid)
  ellipse_y <- b*sin(theta_grid)

  R <- matrix(c(cos(alpha), -sin(alpha), sin(alpha), cos(alpha)), 2, 2)
  val <- cbind(ellipse_x, ellipse_y) %*% R

  val <- t(apply(val, 1, function(x){x + mean_vec}))
}

.ellipse_coverage <- function(dat, multiplier_val, ...){
  stopifnot(ncol(dat) == 2)

  mean_vec <- colMeans(dat)
  cov_mat <- stats::cov(dat)

  e_points <- .ellipse_points(mean_vec, multiplier_val*cov_mat)

  e_hull <- tripack::tri.mesh(epoints[-1,1], epoints[-1,2])
  bool_vec <- tripack::in.convex.hull(e_hull, dat[,1], dat[,2])

  sum(bool_vec)/length(bool_vec)
}

.initial_edges <- function(dat, cluster_labels, target_percentage = 0.8){
  stopifnot(ncol(dat) == 2, all(cluster_labels > 0), all(cluster_labels %% 1 == 0),
            length(cluster_labels) == length(unique(cluster_labels)),
            max(cluster_labels) == length(unique(cluster_labels)))

  # form the ellipses
  k <- max(cluster_labels)
  e_list <- lapply(1:k, function(x){
    dat_subset <- dat[which(cluster_labels == x),]

    scaling_val <- .binary_search(dat_subset, evaluation_func = .ellipse_coverage,
                                  target_val = target_percentage)

    mean_vec <- colMeans(dat)
    cov_mat <- stats::cov(dat)

    .ellipse_points(mean_vec, scaling_val*cov_mat)
  })
}

