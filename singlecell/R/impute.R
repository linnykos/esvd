## Many of the following functions are attributed to scImpute:
##  https://github.com/Vivianstats/scImpute

#' scImpute algorithm
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param drop_idx inidices of \code{dat} to impute
#' @param Kcluster number of clusters
#' @param min_size minimum size of a cluster
#' @param verbose boolean
#' @param weight 0-1 weight factor, where numbers closer to 1 means to rely on newly imputed values
#'
#' @return \code{n} by \code{d} matrix
#' @export
scImpute <- function(dat, drop_idx, Kcluster, min_size = 5, verbose = F,
                      weight = 0.5){
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
                                  weight = weight)
    }
    if(verbose) cat('\n')
  }

  dat2
}


.nnls_impute <- function(cell_vec, neigh_mat, B_vec,
                         max_vec = apply(neigh_mat, 2, max),
                         weight = 1, with_intercept = T){
  stopifnot(weight >= 0 & weight <= 1)
  stopifnot(max(B_vec) <= length(cell_vec))
  stopifnot(!is.matrix(cell_vec))
  stopifnot(length(cell_vec) == ncol(neigh_mat), length(max_vec) == length(cell_vec))
  stopifnot(length(B_vec) == length(unique(B_vec)), length(B_vec) > 1)

  if(length(B_vec) == length(cell_vec)) return(cell_vec)

  # format covariates, also add a constant vector
  if(length(B_vec) <= 15){
    min_val <- min(neigh_mat[neigh_mat > 0])
    x_mat <- t(neigh_mat[,B_vec,drop = F])
    if(with_intercept) x_mat <- cbind(x_mat, min_val)

    nnls_res <- nnls::nnls(A = x_mat, b = cell_vec[B_vec])
    coef_vec <- nnls_res$x

    x_new <- t(neigh_mat[,-B_vec,drop = F])
    if(with_intercept) x_new <- cbind(x_new, min_val)
    y_new <- as.numeric(x_new %*% coef_vec)
  } else {
    x_mat <- t(neigh_mat[,B_vec,drop = F])
    glmnet_res <- glmnet::cv.glmnet(x = x_mat, y = cell_vec[B_vec], family = "gaussian",
                                    intercept=TRUE, nfolds = min(floor(nrow(x_mat)/3), 10))
    x_new <- t(neigh_mat[,-B_vec,drop = F])
    y_new <- as.numeric(stats::predict(glmnet_res, newx = x_new, s = "lambda.min"))
  }

  y_new[y_new < 0] <- 0
  y_new[y_new > max_vec[-B_vec]] <- max_vec[-B_vec][y_new > max_vec[-B_vec]]

  cell_vec2 <- cell_vec
  cell_vec2[-B_vec] <- y_new

  (1-weight)*cell_vec + weight*cell_vec2
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
