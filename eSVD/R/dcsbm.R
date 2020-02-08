## functions to initialize (used from the .project_rank_feasibility function in the initialization.R file)

.spectral_clustering <- function(mat, k){
  svd_res <- .svd_projection(mat, k, factors = T, u_alone = F, v_alone = F)

  u_mat <- t(apply(svd_res$u_mat, 1, function(x){x/.l2norm(x)}))
  v_mat <- t(apply(svd_res$v_mat, 1, function(x){x/.l2norm(x)}))

  row_clustering <- stats::kmeans(u_mat, centers = k, nstart = 20)$cluster
  col_clustering <- stats::kmeans(v_mat, centers = k, nstart = 20)$cluster

  list(row_clustering = row_clustering, col_clustering = col_clustering)
}

# inspired partially by https://arxiv.org/pdf/1612.04717.pdf
.form_prediction_dcsbm <- function(mat, row_clustering, col_clustering){
  n <- nrow(mat); d <- ncol(mat)

  k <- max(row_clustering)

  o_mat <- matrix(NA, k, k)
  for(i in 1:k){
    for(j in 1:k){
      idx_i <- which(row_clustering == i)
      idx_j <- which(col_clustering == j)

      o_mat[i,j] <- sum(mat[idx_i, idx_j])
    }
  }

  o_mat <- o_mat*max(mat)/max(o_mat)

  base_mat <- sapply(1:ncol(mat), function(i){
    o_mat[row_clustering, col_clustering[i]]
  })
  res_mat <- mat/base_mat

  class(res_mat) <- "matrix" #bookkeeping purposes
  nmf_res <- NMF::nmf(res_mat, rank = 1)

  theta_row <- nmf_res@fit@W
  theta_col <- nmf_res@fit@H
  theta_mat <- theta_row %*% theta_col
  svd_res <- .svd_projection(theta_mat, k = 1, factors = T, u_alone = F, v_alone = F)
  theta_row <- abs(svd_res$u_mat)
  theta_col <- abs(svd_res$v_mat)

  theta_row[theta_row >= 1] <- 1
  theta_col[theta_col >= 1] <- 1

  # theta_row <- sapply(1:n, function(i){
  #   sum(mat[i,])/sum(o_mat[row_clustering[i],])
  # })
  # theta_col <- sapply(1:d, function(j){
  #   sum(mat[,j])/sum(o_mat[,col_clustering[j]])
  # })

  list(o_mat = o_mat, theta_row = as.numeric(theta_row), theta_col = as.numeric(theta_col))
}

.dcsbm_projection <- function(mat, k){
  res <- .spectral_clustering(mat, k)
  row_clustering <- res$row_clustering
  col_clustering <- res$col_clustering

  res <- .form_prediction_dcsbm(mat, row_clustering, col_clustering)

  base_mat <- sapply(1:ncol(mat), function(i){
    res$o_mat[row_clustering, col_clustering[i]]
  })

  base_mat <- diag(res$theta_row) %*% base_mat %*% diag(res$theta_col)

  list(mat = base_mat, o_mat = res$o_mat, theta_row = res$theta_row, theta_col = res$theta_col)
}
