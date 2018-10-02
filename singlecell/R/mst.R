.mst <- function(dat, labels){
  stopifnot(nrow(dat) == length(labels))
  stopifnot(all(labels > 0), all(labels %% 1 == 0), length(unique(labels)) == max(labels))

  centers <- t(sapply(1:max(labels), function(x){
    colMeans(dat[which(labels == x),])
  }))

  dist_mat <- as.matrix(stats::dist(centers))
  g <- igraph::graph_from_adjacency_matrix(dist_mat, mode = "undirected",
                                           weighted = TRUE, diag = F)
  list(g = igraph::mst(g), centers = centers)
}

.plot_mst <- function(obj, col_vec, ...){
  stopifnot(nrow(obj$centers) == length(col_vec))

  points(obj$centers[,1], obj$centers[,2], col = col_vec, ...)

  len <- length(col_vec)
  combn_mat <- utils::combn(len, 2)
  for(i in 1:ncol(combn_mat)){
    x <- combn_mat[1,i]; y <- combn_mat[2,i]
    if(igraph::are_adjacent(obj$g, x, y)){
      lines(obj$centers[c(x,y),1], obj$centers[c(x,y),2], col = "black", ...)
    }
  }

  invisible()
}
