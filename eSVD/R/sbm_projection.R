.spectral_clustering <- function(dat, K){
  eigenvectors <- .extract_eigenvectors(dat, K)

  stats::kmeans(eigenvectors, centers=K, nstart=20)$cluster
}

# note: for stability reasons, always ask for 1 more eigenvector than needed
.extract_eigenvectors <- function(dat, K){
  eigen_res <- mgcv::slanczos(dat, K+1)
  idx <- order(abs(eigen_res$values), decreasing = T)[1:K]
  eigen_res$vectors <- eigen_res$vectors[,idx,drop = F]
  eigen_res$values <- eigen_res$values[idx]

  sign_vec <- sign(eigen_res$values)
  eigen_val <- abs(eigen_res$values)

  if(length(eigen_val) == 1) {
    diag_mat <- matrix(sign_vec[1] * sqrt(eigen_val), 1, 1)
  } else {
    diag_mat <- diag(sqrt(eigen_val)) %*% diag(sign_vec)
  }

  eigen_res$vectors %*% diag_mat
}

.form_prediction_sbm <- function(dat, cluster_idx){
  diag(dat) <- NA
  K <- max(cluster_idx)

  b_mat <- matrix(NA, K, K)
  for(i in 1:K){
    for(j in 1:i){
      idx_i <- which(cluster_idx == i)

      if(i != j){
        idx_j <- which(cluster_idx == j)
        if(all(is.na(dat[idx_i, idx_i]))){
          b_mat[i,j] <- 0
        } else {
          b_mat[i,j] <- stats::median(dat[idx_i, idx_j], na.rm = T)
        }
        b_mat[j,i] <- b_mat[i,j]

      } else {
        if(all(is.na(dat[idx_i, idx_i]))){
          b_mat[i,i] <- 0
        } else {
          b_mat[i,i] <- stats::median(dat[idx_i, idx_i], na.rm = T)
        }
      }
    }
  }

  b_mat
}

# optional argument, dat_org which was the original matrix (with missing values)
.sbm_projection <- function(dat, K, dat_org = NA){
  cluster_idx <- .spectral_clustering(dat, K)

  if(all(is.na(dat_org))){
    b_mat <- .form_prediction_sbm(dat, cluster_idx)
  } else {
    b_mat <- .form_prediction_sbm(dat_org, cluster_idx)
  }


  sapply(1:ncol(dat), function(x){
    b_mat[cluster_idx[x], cluster_idx]
  })
}
