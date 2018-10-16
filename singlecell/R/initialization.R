## Many of the following functions are attributed to scImpute:
##  https://github.com/Vivianstats/scImpute

#' Initialization function
#'
#' @param dat matrix
#' @param k numeric
#'
#' @return list
.initialization <- function(dat, k = 2, family = "exponential",
                            max_val = NA){
  stopifnot(length(which(is.na(dat))) == 0)

  idx <- which(dat == 0)
  min_val <- min(dat[which(dat > 0)])
  dat[which(dat == 0)] <- min_val/2
  if(family == "exponential") {
    dat2 <- -1/dat
    direction <- "<="
  } else {
    dat2 <- 1/dat
    direction <- ">="
  }

  svd_res <- svd(dat2)

  if(k == 1) diag_vec <- as.matrix(sqrt(svd_res$d[1:k])) else diag_vec <- sqrt(svd_res$d[1:k])
  u_mat <- svd_res$u[,1:k,drop = F] %*% diag(diag_vec)
  v_mat <- svd_res$v[,1:k,drop = F] %*% diag(diag_vec)

  # project v back into positive space based on u
  for(j in 1:nrow(v_mat)){
    v_mat[j,] <- .projection_l1(v_mat[j,], u_mat, which(!is.na(dat[,j])),
                                direction = direction, other_bound = max_val)
  }

  pred_mat <- u_mat %*% t(v_mat)
  if(family == "exponential") {
    stopifnot(all(pred_mat[which(!is.na(dat))] < 0))
  } else {
    stopifnot(all(pred_mat[which(!is.na(dat))] > 0))
  }

  list(u_mat = u_mat, v_mat = v_mat)
}

.find_neighbors_impute <- function(dat, var_thres = 0.6, Kcluster = 2){
  stopifnot(length(which(is.na(dat))) == 0)

  pca_res <- stats::prcomp(dat)
  eigs <- (pca_res$sdev)^2
  var_cum <- cumsum(eigs)/sum(eigs)
  if(max(var_cum) <= var_thres){
    npc <- length(var_cum)
  } else {
    npc <- which.max(var_cum > var_thres)
  }

  if (npc < 3){ npc = 3 }
  mat_pcs <- t(pca_res$x[, 1:npc])

  J <- ncol(mat_pcs) # number of cells
  dist_cells_list <- lapply(1:J, function(id1){
    d <- sapply(1:id1, function(id2){
      sse <- sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
      sqrt(sse)
    })
    c(d, rep(0, J-id1))
  })
  dist_cells <- matrix(0, nrow = J, ncol = J)
  for(cellid in 1:J){
    dist_cells[cellid, ] = dist_cells_list[[cellid]]
  }
  dist_cells <- dist_cells + t(dist_cells)

  non_out <- 1:J
  spec_res <- kernlab::specc(t(mat_pcs[, non_out]),
                             centers = Kcluster,
                             kernel = "rbfdot")
  nbs <- rep(NA, J)
  nbs[non_out] <- spec_res

  nbs
}

.nnls_impute <- function(cell_vec, neigh_mat, B_vec,
                         max_vec = apply(neigh_mat, 2, max),
                         max_time = 60){
  stopifnot(max(B_vec) <= length(cell_vec))
  stopifnot(!is.matrix(cell_vec))
  stopifnot(length(cell_vec) == ncol(neigh_mat), length(max_vec) == length(cell_vec))
  stopifnot(length(B_vec) == length(unique(B_vec)), length(B_vec) > 1)

  if(length(B_vec) == length(cell_vec)) return(cell_vec)

  # format covariates, also add a constant vector
  min_val <- min(cell_vec[cell_vec > 0])
  x_mat <- t(neigh_mat[,B_vec,drop = F])
  x_mat <- cbind(x_mat, min_val)

  nnls_res <- nnls::nnls(A = x_mat, b = cell_vec[B_vec])
  coef_vec <- nnls_res$x
  y_new <- as.numeric(cbind(t(neigh_mat[,-B_vec,drop = F]), min_val) %*% coef_vec)

  y_new[y_new > max_vec[-B_vec]] <- max_vec[-B_vec][y_new > max_vec[-B_vec]]

  cell_vec2 <- cell_vec
  cell_vec2[-B_vec] <- y_new

  cell_vec2
}

.scImpute <- function(dat, drop_idx, Kcluster, min_size = 5, max_time = 60, verbose = F){
  stopifnot(length(which(is.na(dat))) == 0)
  if(length(drop_idx) == 0) return(dat)

  neigh_vec <- .find_neighbors_impute(dat, Kcluster = Kcluster)
  neigh_list <- lapply(1:Kcluster, function(k){which(neigh_vec == k)})
  stopifnot(all(sapply(neigh_list, length) >= min_size))

  dat2 <- dat
  dat2[drop_idx] <- NA

  for(k in neigh_list){
    for(i in 1:length(k)){
      if(verbose && i %% floor(length(k)/10) == 0) cat('*')
      keep_idx <- which(!is.na(dat2[k[i],]))
      dat2[k[i],] <- .nnls_impute(dat[k[i],], dat[setdiff(k, k[i]),,drop = F], keep_idx,
                               max_time = max_time)
    }
    if(verbose) cat('\n')
  }

  dat2
}
