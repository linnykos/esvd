slingshot_3dplot <- function(dat, cluster_labels, bg_col_vec,
                             cluster_col_vec, breaks,
                             bg_cex = 0.4, cluster_cex = 2,
                             curves = NA, ...){
  stopifnot(length(bg_col_vec) == length(breaks) - 1 &
              length(cluster_col_vec) == length(breaks) - 1)

  cluster_center <- .compute_cluster_center(dat, .construct_cluster_matrix(cluster_labels))

  plot3D::scatter3D(x = dat[,1], y = dat[,2], z = dat[,3],
                    surface = FALSE, colvar = cluster_labels,
                    cex = bg_cex,
                    breaks = breaks, col = bg_col_vec, colkey = F, ...)
  plot3D::points3D(cluster_center[,1], cluster_center[,2], cluster_center[,3],
                   colvar = 1:max(cluster_labels),
                   breaks = breaks, add = T, col = cluster_col_vec,
                   cex = cluster_cex, colkey = F, ...)

  if(!all(is.na(curves))){
    for(k in 1:length(curves)){
      ord <- curves[[k]]$ord
      plot3D::lines3D(x = curves[[k]]$s[ord, 1],
                      y = curves[[k]]$s[ord, 2],
                      z = curves[[k]]$s[ord, 3],
                      add = T, colkey = F, col = "black", ...)
    }
  }

  invisible()
}

#' Construct the 3D tube for each curve
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param radius a numeric
#' @param len a numeric
#'
#' @return a list
#' @export
construct_3d_tube <- function(dat, radius, len = 20){
  stopifnot(nrow(dat) > 2)

  n <- nrow(dat)
  circle_list <- lapply(1:nrow(dat), function(x){
    direction <- .find_adjacent_directions(dat, x)
    res <- .find_basis_vectors(direction)
    basis_vec1 <- res$vec1; basis_vec2 <- res$vec2

    .construct_3d_circle(dat[x,], radius, basis_vec1, basis_vec2, len = len)
  })

  # apply the correction
  for(i in 2:length(circle_list)){
    circle_list[[i]] <- .find_correct_orientation(circle_list[[i-1]], circle_list[[i]])
  }

  # form the matrices
  x_mat <- matrix(NA, nrow = nrow(circle_list[[1]]), ncol = length(circle_list))
  y_mat <- x_mat; z_mat <- x_mat
  for(i in 1:length(circle_list)){
    x_mat[,i] <- circle_list[[i]][,1]
    y_mat[,i] <- circle_list[[i]][,2]
    z_mat[,i] <- circle_list[[i]][,3]
  }

  list(x_mat = x_mat, y_mat = y_mat, z_mat = z_mat)
}


.find_adjacent_directions <- function(dat, idx){
  stopifnot(idx >= 1)
  stopifnot(idx <= nrow(dat))

  idx_vec <- sort(unique(pmax(pmin(idx + c(-1,0,1), nrow(dat)), 1)))
  stopifnot(length(idx_vec) %in% c(2,3))

  res <- rowMeans(sapply(1:(length(idx_vec)-1), function(i){
    vec <- dat[idx_vec[i],] - dat[idx_vec[i+1],]
    vec/.l2norm(vec)
  }))

  res/.l2norm(res)
}

.projection_matrix <- function(vec, mat){
  stopifnot(length(vec) == nrow(mat), nrow(mat) >= ncol(mat))
  n <- length(vec)
  proj_mat <- diag(n) - mat %*% solve(t(mat) %*% mat) %*% t(mat)

  as.numeric(proj_mat %*% vec)
}

.find_basis_vectors <- function(direction){
  stopifnot(length(direction) == 3)

  mat <- cbind(direction, matrix(stats::rnorm(6), ncol = 2, nrow = 3))
  for(i in 2:3) {
    mat[,i] <- .projection_matrix(mat[,i], mat[,1:(i-1), drop = F])
    mat[,i] <- mat[,i]/.l2norm(mat[,i])
  }

  list(vec1 = mat[,2], vec2 = mat[,3])
}

.construct_3d_circle <- function(dat_vec, radius, basis_vec1, basis_vec2,
                                 len = 20){
  seq_vec <- seq(0, 2*pi, length.out = len)
  x <- radius*cos(seq_vec)
  y <- radius*sin(seq_vec)

  mat <- t(cbind(basis_vec1, basis_vec2) %*% rbind(x,y))
  for(i in 1:nrow(mat)){
    mat[i,] <- mat[i,]+dat_vec
  }

  mat
}

.construct_all_circles <- function(dat, radius, len = 20){
  n <- nrow(dat)
  lapply(1:n, function(x){
    direction <- .find_adjacent_directions(dat, x)
    res <- .find_basis_vectors(direction)
    basis_vec1 <- res$vec1; basis_vec2 <- res$vec2
    .construct_3d_circle(dat[x,], radius, basis_vec1, basis_vec2, len = len)
  })
}

.find_correct_orientation <- function(mat1, mat2){
  stopifnot(ncol(mat1) == 3, ncol(mat2) == 3, nrow(mat1) == nrow(mat2))
  n <- nrow(mat1)

  func <- function(vec, mat){
    dist_vec <- sapply(1:nrow(mat), function(x){
      .l2norm(vec - mat[x,])
    })

    c(min(dist_vec), which.min(dist_vec))
  }

  res <- sapply(1:nrow(mat1), function(x){
    func(mat1[x,], mat2)
  })

  anchor1 <- which.min(res[1,])
  anchor2 <- res[2,anchor1]

  direction_forward <- (length(which(diff(res[2,]) > 0)) > length(which(diff(res[2,]) < 0)))

  # forward direction
  indices_list <- .subset_indices(n, anchor1, anchor2, direction_forward)
  mat2[indices_list$to,] <- mat2[indices_list$from,]

  mat2
}

.subset_indices <- function(n, anchor1, anchor2, direction_forward = T){
  vec1 <- anchor1:(anchor1+n-1)
  if(direction_forward){
    vec2 <- anchor2:(anchor2+n-1)
  } else {
    vec2 <- anchor2:(anchor2-n+1)
  }

  stopifnot(length(vec1) == n, length(vec2) ==n)

  mat <- rbind(vec1, vec2)
  mat[mat > n] <- mat[mat > n] - n
  mat[mat < 1] <- mat[mat < 1] + n

  stopifnot(all(mat >= 1), all(mat <= n))

  list(to = mat[1,], from = mat[2,])
}







