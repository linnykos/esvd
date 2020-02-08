## functions to initialize (used from the .project_rank_feasibility function in the initialization.R file)

.spectral_clustering <- function(mat, k){
  svd_res <- .svd_projection(mat, k, factors = T, u_alone = F, v_alone = F)

  u_mat <- t(apply(svd_res$u_mat, 1, function(x){x/.l2norm(x)}))
  v_mat <- t(apply(svd_res$v_mat, 1, function(x){x/.l2norm(x)}))

  # make sure there are more than k distinct points
  if(length(unique(u_mat[,1])) <= k) u_mat[,1] <- u_mat[,1] + rnorm(nrow(u_mat), sd = 0.01)
  if(length(unique(v_mat[,1])) <= k) v_mat[,1] <- v_mat[,1] + rnorm(nrow(v_mat), sd = 0.01)

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

  # preparing to estimate the node-specific factors
  o_mat <- o_mat*max(mat)/max(o_mat)

  base_mat <- sapply(1:ncol(mat), function(i){
    o_mat[row_clustering, col_clustering[i]]
  })
  res_mat <- mat/base_mat

  # fitting the NMF
  class(res_mat) <- "matrix" #bookkeeping purposes
  nmf_res <- NMF::nmf(res_mat, rank = 1)

  # adjusting the NMF fit
  theta_row <- nmf_res@fit@W
  theta_col <- nmf_res@fit@H
  theta_mat <- theta_row %*% theta_col
  svd_res <- .svd_projection(theta_mat, k = 1, factors = T, u_alone = F, v_alone = F)
  theta_row <- abs(svd_res$u_mat)
  theta_col <- abs(svd_res$v_mat)

  # thresholding the theta's to ensure the predictions are within the range given by mat
  theta_row[theta_row >= 1] <- 1
  theta_col[theta_col >= 1] <- 1

  list(o_mat = o_mat, theta_row = as.numeric(theta_row), theta_col = as.numeric(theta_col))
}

.dcsbm_projection <- function(mat, k){
  res <- .spectral_clustering(mat, k)
  row_clustering <- res$row_clustering
  col_clustering <- res$col_clustering

  sign_val <- ifelse(all(mat <= 0), -1, 1)
  mat <- sign_val * mat
  res <- .form_prediction_dcsbm(mat, row_clustering, col_clustering)

  base_mat <- sapply(1:ncol(mat), function(i){
    res$o_mat[row_clustering, col_clustering[i]]
  })

  base_mat <- diag(res$theta_row) %*% base_mat %*% diag(res$theta_col)

  list(mat = sign_val*base_mat, o_mat = res$o_mat, theta_row = res$theta_row, theta_col = res$theta_col)
}
