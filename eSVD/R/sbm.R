## functions to initialize (used from the .project_rank_feasibility function in the initialization.R file)

#' Spectral clustering
#'
#' @param dat matrix where spectral clustering will be applied to both the rows and columns
#' @param k integer
#'
#' @return a list \code{row_clustering} and \code{col_clustering} for the respective clusterings
.spectral_clustering <- function(dat, k){
  stopifnot(k <= min(dim(dat)))

  svd_res <- .svd_projection(dat, k, factors = T, u_alone = F, v_alone = F)

  row_clustering <- stats::kmeans(svd_res$u_mat, centers = k, nstart = 20)$cluster
  col_clustering <- stats::kmeans(svd_res$v_mat, centers = k, nstart = 20)$cluster

  list(row_clustering = row_clustering, col_clustering = col_clustering)
}

#' Form biclustering predictions based on clustering of rows and columns
#'
#' It is assumed that both \code{row_clustering} and \code{col_clustering} have the
#' same number of clusters
#'
#' @param dat a \code{n} by \code{p} matrix
#' @param row_clustering a clustering of \code{n} elements
#' @param col_clustering a clustering of \code{p} elements
#'
#' @return a \code{k} by \code{k} matrix
.form_prediction_sbm <- function(dat, row_clustering, col_clustering){
  stopifnot(all(row_clustering > 0), all(row_clustering %% 1 == 0), length(row_clustering) == nrow(dat),
            all(col_clustering > 0), all(col_clustering %% 1 == 0), length(col_clustering) == ncol(dat),
            max(row_clustering) == max(col_clustering))

  k <- max(row_clustering)

  b_mat <- matrix(NA, k, k)
  for(i in 1:k){
    for(j in 1:k){
      idx_i <- which(row_clustering == i)
      idx_j <- which(col_clustering == j)

      b_mat[i,j] <- stats::median(dat[idx_i, idx_j], na.rm = T)
    }
  }

  b_mat
}

#' Projection based on the stochastic block model
#'
#' @param dat matrix
#' @param k integer
#'
#' @return matrix of the same dimensions as \code{dat}, but has only at most \code{k^2}
#' unique values and is guaranteed to be rank \code{k}
.sbm_projection <- function(dat, k){
  res <- .spectral_clustering(dat, k)

  b_mat <- .form_prediction_sbm(dat, res$row_clustering, res$col_clustering)

  sapply(1:ncol(dat), function(i){
    b_mat[res$row_clustering, res$col_clustering[i]]
  })
}
