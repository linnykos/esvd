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
.get_lineages <- function(dat, cluster_labels, starting_cluster, knn = NA,
                          remove_outlier = T, percentage = 0.05){
  ### formatting
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  stopifnot(is.matrix(dat), nrow(dat) == nrow(cluster_mat))

  ### get the connectivity matrix
  centers <- .compute_cluster_center(dat, cluster_mat)
  dat_augment <- rbind(centers, dat)

  ### construct the k-nearest neighbor graph
  knn_graph <- .determine_knn(dat_augment, knn)

  if(remove_outlier){
    idx <- .remove_outliers(knn_graph, cluster_labels, percentage = percentage)
    knn_graph <- .determine_knn(dat_augment[c(1:k, idx+k),], knn)
  }

  ### construct the spt
  spt_graph <- .construct_spt(knn_graph, k = k, starting_cluster = starting_cluster)

  ### identify lineages (paths through trees)
  lineages <- .construct_lineages(spt_graph, starting_cluster = starting_cluster)

  lineages
}

#############

#' Construst cluster matrix from cluster labels
#'
#' @param cluster_labels  vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels, na.rm = T)}. Can include \code{NA}
#'
#' @return A 0-1 matrix with \code{length(cluster_labels)} rows
#' and \code{max(cluster_labels)} columns
.construct_cluster_matrix <- function(cluster_labels){
  idx <- !is.na(cluster_labels)
  stopifnot(all(cluster_labels[idx] > 0), all(cluster_labels[idx] %% 1 == 0))
  stopifnot(length(unique(cluster_labels[idx])) == max(cluster_labels[idx]))

  k <- max(cluster_labels, na.rm = T)
  n <- length(cluster_labels)

  mat <- sapply(1:k, function(x){
    tmp <- rep(0, n)
    tmp[which(cluster_labels == x)] <- 1
    tmp
  })

  colnames(mat) <- seq_len(k)
  mat
}

#' Compute the cluster centers
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param cluster_mat a 0-1 matrix that is \code{n} by \code{k}
#'
#' @return a \code{k} by \code{d} matrix
.compute_cluster_center <- function(dat, cluster_mat){
  mat <- t(sapply(1:ncol(cluster_mat), function(x){
    idx <- which(cluster_mat[,x] == 1)
    colMeans(dat[idx,,drop=F])
  }))
  rownames(mat) <- colnames(cluster_mat)
  mat
}

#' Construct the K-nearest neighbor graph based on Euclidean distance
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param knn positive integer
#'
#' @return \code{igraph} graph object
.construct_knn_graph <- function(dat, knn = 5){
  n <- nrow(dat)
  dist_mat <- as.matrix(stats::dist(dat, method = "euclidean"))
  adj_mat <- sapply(1:n, function(x){
    vec <- dist_mat[,x]; vec[x] <- Inf
    idx <- order(vec, decreasing = F)[1:knn]
    adj_vec <- rep(0, n)
    adj_vec[idx] <- 1
    adj_vec
  })

  adj_mat <- ceiling((adj_mat + t(adj_mat))/2)
  igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
}

#' Construct the shortest path tree graph from KNN graph
#'
#' This is currently hard-coded to work for only one \code{starting_cluster}.
#'
#' Note: The squaring of \code{dist_mat} is arbitrary, currently used to
#' encourage paths through other clusters.
#'
#' @param knn_graph \code{igraph} object
#' @param k positive integer for number of clusters
#' @param starting_cluster the "origin" cluster that all the lineages will start
#' from
#'
#' @return \code{igraph} object representing the shortest path tree
.construct_spt <- function(knn_graph, k, starting_cluster){
  stopifnot(starting_cluster <= k)

  # construct distance graph
  dist_mat <- igraph::distances(knn_graph, v = 1:k, to = 1:k, mode = "all")
  dist_mat[is.infinite(dist_mat)] <- 0
  dist_mat <- dist_mat^2
  dist_graph <- igraph::graph_from_adjacency_matrix(dist_mat, weighted = T,
                                                    mode = "undirected")

  # enumerate shortest paths and then make it into a graph
  lis <- igraph::shortest_paths(dist_graph, from = starting_cluster, output = "vpath")
  adj_mat <- matrix(0, k, k)
  for(vpath in lis$vpath){
    if(length(vpath) == 1) next()
    for(i in 2:length(vpath)){
      adj_mat[vpath[i-1], vpath[i]] <- 1
    }
  }
  adj_mat <- ceiling((adj_mat + t(adj_mat))/2)
  igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
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
.construct_lineages <- function(spt_graph, starting_cluster){
  # find all leaf nodes (aside from starting_cluster) and trace its path
  deg_vec <- igraph::degree(spt_graph)
  leaf_idx <- which(deg_vec == 1)
  leaf_idx <- leaf_idx[!leaf_idx %in% starting_cluster]

  stopifnot(length(leaf_idx) > 0)
  lineages <- igraph::shortest_paths(spt_graph, from = starting_cluster,
                                     to = leaf_idx, output = "vpath")$vpath
  for(i in 1:length(lineages)){
    lineages[[i]] <- as.numeric(lineages[[i]])
  }

  names(lineages) <- paste('Lineage',seq_along(lineages),sep='')

  lineages
}

#' Remove outliers based on shortest path centrality
#'
#' This function relies on the specific format of \code{graph}.
#' Specifically, the first \code{max(cluster_labels, na.rm=T)} nodes
#' in \code{graph} correspond to the clutser centers and are not
#' reflected with \code{cluster_labels}
#'
#' @param graph \code{igraph} object
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels, na.rm = T)}. Can include \code{NA}
#' @param percentage percentage of vertices to remove that have
#' a cluster label
#'
#' @return indices (of the original dataset, reflected by \code{cluster_labels})
#' that should be kept
.remove_outliers <- function(graph, cluster_labels, percentage = 0.05){
  stopifnot(class(graph) == "igraph")
  centrality <- igraph::betweenness(graph)

  k <- max(cluster_labels, na.rm = T)

  # remove outliers based on IQR
  iqr_mat <- do.call(rbind, lapply(1:k, function(x){
    idx <- which(cluster_labels == x)
    iqr_val <- stats::IQR(centrality[idx+k])
    med_val <- stats::median(centrality[idx+k])

    cbind(idx, x, (centrality[idx+k]-med_val)/iqr_val)
  }))
  iqr_mat <- iqr_mat[order(iqr_mat[,ncol(iqr_mat)], decreasing = T),]
  n <- nrow(iqr_mat)
  remove_percentage <- round(percentage * n)

  iqr_mat[-c(1:remove_percentage),1]
}

.determine_knn <- function(dat_augment, knn){
  if(is.na(knn)){
    knn <- 1
    while(TRUE){
      knn_graph <- .construct_knn_graph(dat_augment, knn = knn)
      if(igraph::components(knn_graph)$no == 1) break()
      knn <- knn + 1
    }
  } else {
    knn_graph <- .construct_knn_graph(dat_augment, knn = knn)
    stopifnot(igraph::components(knn_graph)$no == 1)
  }

  knn_graph
}
