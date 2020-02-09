## functions to initialize (used from the .project_rank_feasibility function in the initialization.R file)

.spectral_clustering <- function(dat, k){
  svd_res <- .svd_projection(dat, k, factors = T, u_alone = F, v_alone = F)

  row_clustering <- stats::kmeans(svd_res$u_mat, centers = k, nstart = 20)$cluster
  col_clustering <- stats::kmeans(svd_res$v_mat, centers = k, nstart = 20)$cluster

  list(row_clustering = row_clustering, col_clustering = col_clustering)
}

.form_prediction_sbm <- function(dat, row_clustering, col_clustering){
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

.sbm_projection <- function(dat, k){
  res <- .spectral_clustering(dat, k)

  b_mat <- .form_prediction_sbm(dat, res$row_clustering, res$col_clustering)

  sapply(1:ncol(dat), function(i){
    b_mat[res$row_clustering, res$col_clustering[i]]
  })
}
