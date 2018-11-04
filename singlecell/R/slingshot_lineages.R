# code adapted from https://github.com/kstreet13/slingshot

#' Title
#'
#' Note: I removed the functionality to explicitly label a starting or
#' ending cluster (might put back in later?).
#'
#' Note: I removed the Omega parameter, which, to my understanding,
#' controls if a lineage (i.e. tree) is split into two separate lineages.
#'
#' @param dat
#' @param cluster_labels
#' @param knn
#'
#' @return
#' @export
#'
#' @examples
.get_lineages <- function(dat, cluster_labels, knn = 5){
  ### formatting
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  stopifnot(is.matrix(dat), nrow(dat) == nrow(cluster_mat))

  ### get the connectivity matrix
  centers <- .compute_cluster_center(dat, cluster_mat)
  dat_augment <- rbind(centers, dat)

  ### construct the k-nearest neighbor graph
  knn_graph <- .construct_knn_graph(dat_augment, knn = knn)

  # ### get pairwise cluster distance matrix
  # D <- .pairwise_cluster_distance(dat, cluster_mat)
  #
  # ### set omega
  # if(is.na(omega)) omega <- max(D) + 1
  # D <- rbind(D, rep(omega, ncol(D)))
  # D <- cbind(D, c(rep(omega, ncol(D)), 0))

  ### draw MST on cluster centers + OMEGA
  forest <- .construct_mst_slingshot(D, cluster_mat, end_clus)

  ### identify sub-trees
  trees <- .identify_trees(forest)
  ntree <- length(trees)

  ### identify lineages (paths through trees)
  lineages <- .construct_lineages(trees, forest, start_clus)
  lineageControl <- .format_lineagecontrol(lineages, D, cluster_mat, start_clus, end_clus)

  list(dat = dat, cluster_mat = cluster_mat,
       lineages = lineages, adjacency = forest,
       slingParams = lineageControl)
}

#############

#' Construst cluster matrix from cluster labels
#'
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}
#'
#' @return A 0-1 matrix with \code{length(cluster_labels)} rows
#' and \code{max(cluster_labels)} columns
.construct_cluster_matrix <- function(cluster_labels){
  stopifnot(all(cluster_labels > 0), all(cluster_labels %% 1 == 0))
  stopifnot(length(unique(cluster_labels)) == max(cluster_labels))

  k <- max(cluster_labels)
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
  t(sapply(1:ncol(cluster_mat), function(x){
    idx <- which(cluster_mat[,x] == 1)
    colMeans(dat[idx,,drop=F])
  }))
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


.construct_mst_slingshot <- function(D, cluster_mat, end_clus = NA){
  clusters <- 1:ncol(cluster_mat)
  nclus <- ncol(cluster_mat)

  if(!any(is.na(end_clus))){
    end_idx <- which(clusters %in% end_clus)
    mstree <- ape::mst(D[-end_clus, -end_clus, drop = FALSE])
  }else{
    mstree <- ape::mst(D)
  }

  # (add in endpoint clusters)
  if(!any(is.na(end_clus))){
    forest <- D
    forest[forest != 0] <- 0
    forest[-end.idx, -end.idx] <- mstree
    for(cluster_id in end.clus){
      cl_idx <- which(clusters == cluster_id)
      dists <- D[! rownames(D) %in% end_clus, cl_idx]
      # get closest non-endpoint cluster
      closest <- names(dists)[which.min(dists)]
      closest_idx <- which.max(clusters == closest)
      forest[cl_idx, closest_idx] <- 1
      forest[closest_idx, cl_idx] <- 1
    }
  }else{
    forest <- mstree
  }
  # remove OMEGA
  forest <- forest[seq_len(nclus), seq_len(nclus), drop = FALSE]
  rownames(forest) <- clusters
  colnames(forest) <- clusters

  forest
}

.identify_trees <- function(forest){
  subtrees <- subtrees.update <- forest
  diag(subtrees) <- 1
  while(sum(subtrees.update) > 0){
    subtrees.new <- apply(subtrees,2,function(col){
      rowSums(subtrees[,as.logical(col), drop=FALSE]) > 0
    })
    subtrees.update <- subtrees.new - subtrees
    subtrees <- subtrees.new
  }
  subtrees <- unique(subtrees)
  trees <- lapply(seq_len(nrow(subtrees)),function(ri){
    colnames(forest)[subtrees[ri,]]
  })

  trees[order(vapply(trees,length,0),decreasing = TRUE)]
}

.construct_lineages <- function(trees, forest, start_clus = NA){
  lineages <- list()

  for(tree in trees){
    if(length(tree) == 1){
      lineages[[length(lineages)+1]] <- tree
      next
    }

    tree_ind <- rownames(forest) %in% tree
    tree_graph <- forest[tree_ind, tree_ind, drop = FALSE]
    degree <- rowSums(tree_graph)
    g <- igraph::graph.adjacency(tree_graph, mode="undirected")

    # if you have starting cluster(s) in this tree, draw lineages
    # to each leaf
    if(!is.na(any(start_clus)) && sum(start_clus %in% tree) > 0){
      starts <- start_clus[start_clus %in% tree]
      ends <- rownames(tree_graph)[
        degree == 1 & ! rownames(tree_graph) %in% starts]
      for(st in starts){
        paths <- igraph::shortest_paths(g, from = st, to = ends,
                                mode = 'out',
                                output = 'vpath')$vpath
        for(p in paths){
          lineages[[length(lineages)+1]] <- names(p)
        }
      }
    } else {
      # else, need a criteria for picking root
      # highest average length (~parsimony)
      leaves <- rownames(tree_graph)[degree == 1]
      avg_lineage_length <- vapply(leaves, function(l){
        ends <- leaves[leaves != l]
        paths <- igraph::shortest_paths(g, from = l, to = ends,
                                        mode = 'out',
                                        output = 'vpath')$vpath
        mean(vapply(paths, length, 0))
      }, 0)
      st <- names(avg_lineage_length)[
        which.max(avg_lineage_length)]
      ends <- leaves[leaves != st]
      paths <- igraph::shortest_paths(g, from = st, to = ends,
                                      mode = 'out',
                                      output = 'vpath')$vpath
      for(p in paths){
        lineages[[length(lineages)+1]] <- names(p)
      }
    }
  }

  lineages <- lineages[order(vapply(lineages, length, 0),
                             decreasing = TRUE)]
  names(lineages) <- paste('Lineage',seq_along(lineages),sep='')

  lineages
}

.format_lineagecontrol <- function(lineages, D, cluster_mat, start_clus, end_clus){
  nclus <- ncol(cluster_mat)

  lineageControl <- list()
  first <- unique(vapply(lineages,function(l){ l[1] },''))
  last <- unique(vapply(lineages,function(l){ l[length(l)] },''))

  lineageControl$start_clus <- first
  lineageControl$end_clus <- last

  start_given <- first %in% start_clus
  end_given <- last %in% end_clus
  lineageControl$start_given <- start_given
  lineageControl$end_given <- end_given

  lineageControl$dist <- D[seq_len(nclus),seq_len(nclus),
                           drop = FALSE]

  lineageControl
}
