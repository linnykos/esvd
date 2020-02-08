## functions to initialize (used from the .project_rank_feasibility function in the initialization.R file)

.spectral_clustering <- function(dat, k){
  svd_res <- .svd_projection(dat, k, factors = T, u_alone = F, v_alone = F)

  u_mat <- t(apply(svd_res$u_mat, 1, function(x){x/.l2norm(x)}))
  v_mat <- t(apply(svd_res$v_mat, 1, function(x){x/.l2norm(x)}))

  row_clustering <- stats::kmeans(u_mat, centers = k, nstart = 20)$cluster
  col_clustering <- stats::kmeans(v_mat, centers = k, nstart = 20)$cluster

  list(row_clustering = row_clustering, col_clustering = col_clustering)
}

# from https://arxiv.org/pdf/1612.04717.pdf
.form_prediction_dcsbm <- function(dat, row_clustering, col_clustering){
  n <- nrow(dat); d <- ncol(dat)

  k <- max(row_clustering)

  o_mat <- matrix(NA, k, k)
  for(i in 1:k){
    for(j in 1:k){
      idx_i <- which(row_clustering == i)
      idx_j <- which(col_clustering == j)

      o_mat[i,j] <- sum(dat[idx_i, idx_j], na.rm = T)
    }
  }

  theta_row <- sapply(1:n, function(i){
    sum(dat[i,])/sum(o_mat[row_clustering[i],])
  })
  theta_col <- sapply(1:d, function(j){
    sum(dat[,j])/sum(o_mat[,col_clustering[j]])
  })

  list(o_mat = o_mat, theta_row = theta_row, theta_col = theta_col)
}

.dcsbm_projection <- function(dat, k){
  res <- .spectral_clustering(dat, k)
  row_clustering <- res$row_clustering
  col_clustering <- res$col_clustering

  res <- .form_prediction_dcsbm(dat, row_clustering, col_clustering)

  base_mat <- sapply(1:ncol(dat), function(i){
    res$o_mat[row_clustering, col_clustering[i]]
  })

  base_mat <- diag(res$theta_row) %*% base_mat %*% diag(res$theta_col)

  list(mat = base_mat, o_mat = res$o_mat, theta_row = res$theta_row, theta_col = res$theta_col)
}
