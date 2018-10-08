# code adapted from https://github.com/kstreet13/slingshot

.get_lineages <- function(dat, cluster_labels, start_clus = NA, end_clus = NA,
                          omega = NA){
  ### formatting
  cluster_mat <- .format_clusterlabels(cluster_labels)
  stopifnot(is.matrix(dat), nrow(dat) == nrow(cluster_mat))
  if(!any(is.na(start_clus))){
    start_clus <- as.character(start_clus)
  }
  if(!any(is.na(end_clus))){
    end_clus <- as.character(end_clus)
  }

  ### get the connectivity matrix
  centers <- .compute_clustercenter(dat, cluster_mat)

  ### get pairwise cluster distance matrix
  D <- .pairwise_cluster_distance(dat, cluster_mat)

  ### set omega
  if(is.na(omega)) omega <- max(D) + 1
  D <- rbind(D, rep(omega, ncol(D)))
  D <- cbind(D, c(rep(omega, ncol(D)), 0))

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

.format_clusterlabels <- function(cluster_labels){
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

.compute_clustercenter <- function(dat, cluster_mat){
  t(vapply(1:ncol(cluster_mat), function(cluster_id){
    w <- cluster_mat[,cluster_id]
    matrixStats::colWeightedMeans(dat, w = w)
  }, rep(0,ncol(dat))))
}

.dist_clusters_full <- function(dat, w1, w2){
  if(length(w1) != nrow(dat) | length(w2) != nrow(dat)){
    stop("Reduced dimensional matrix and weights vector contain different
         numbers of points.")
  }

  mu1 <- matrixStats::colWeightedMeans(dat, w = w1)
  mu2 <- matrixStats::colWeightedMeans(dat, w = w2)
  diff <- mu1 - mu2
  s1 <- stats::cov.wt(dat, wt = w1)$cov
  s2 <- stats::cov.wt(dat, wt = w2)$cov

  as.numeric(t(diff) %*% solve(s1 + s2) %*% diff)
}

.pairwise_cluster_distance <- function(dat, cluster_mat){
  clusters <- 1:ncol(cluster_mat)
  nclus <- ncol(cluster_mat)

  D <- as.matrix(vapply(clusters, function(cluster_id1){
    vapply(clusters, function(cluster_id2){
      w1 <- cluster_mat[,cluster_id1]
      w2 <- cluster_mat[,cluster_id2]

      .dist_clusters_full(dat, w1, w2)
    },0)
  }, rep(0, nclus)))

  rownames(D) <- as.character(clusters)
  colnames(D) <- as.character(clusters)

  D
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
