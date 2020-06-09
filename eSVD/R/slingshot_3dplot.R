#' 3D plotting function
#'
#' Relies heavily on the \code{plot3D} package
#'
#' @param dat \code{n} by \code{k} (required to be 3 or larger) dataset, of which the first 3 columns
#' are plotted
#' @param cluster_labels numeric of the same length as \code{nrow(dat)}, where the values are consecutive
#' integers from 1 to \code{max(cluster_labels)}
#' @param bg_col_vec vector of colors of the same length as \code{max(cluster_labels)}
#' @param bg_cex numeric for plotting, denoting the size of each point (i.e., row) in \code{dat}
#' @param cluster_center location of the cluster centers, possibly \code{NA}. Otherwise, it is a
#' numeric matrix of \code{max(cluster_labels)} rows and at least 3 columns.
#' @param center_cex numeric for plotting, denoting the size of each cluster center in \code{cluster_center}
#' @param center_col_vec  vector of colors of the same length as \code{max(cluster_labels)}
#' @param ... additional plotting parameter
#'
#' @return no return. A plot is produced
#' @export
slingshot_3dplot <- function(dat, cluster_labels, bg_col_vec, bg_cex = 0.4,
                             cluster_center = NA, center_cex = 2, center_col_vec = NA, ...){
  stopifnot(length(bg_col_vec) == max(cluster_labels), ncol(dat) >= 3,
            length(cluster_labels) == nrow(dat), is.numeric(cluster_labels))
  stopifnot(all(cluster_labels >= 1), all(cluster_labels %% 1 == 0),
            length(unique(cluster_labels)) == max(cluster_labels))
  if(!all(is.na(cluster_center))){
    stopifnot(ncol(cluster_center) >= 3, nrow(cluster_center) == max(cluster_labels),
              length(center_col_vec) == max(cluster_labels))

  }
  breaks <- seq(0.5, max(cluster_labels)+1, by = 1)

  plot3D::scatter3D(x = dat[,1], y = dat[,2], z = dat[,3],
                    surface = FALSE, colvar = cluster_labels,
                    cex = bg_cex,
                    breaks = breaks, col = bg_col_vec, colkey = F, ...)

  if(!all(is.na(cluster_center))){
    center_labels <- 1:nrow(cluster_center)

    plot3D::points3D(cluster_center[,1], cluster_center[,2], cluster_center[,3],
                     colvar = center_labels,
                     breaks = breaks, add = T, col = rep("black", 13),
                     cex = 1.5*center_cex, colkey = F, ...)
    plot3D::points3D(cluster_center[,1], cluster_center[,2], cluster_center[,3],
                     colvar = center_labels,
                     breaks = breaks, add = T, col = center_col_vec,
                     cex = center_cex, colkey = F, ...)
    plot3D::points3D(cluster_center[,1], cluster_center[,2], cluster_center[,3],
                     colvar = center_labels,
                     breaks = breaks, add = T, col = center_col_vec,
                     cex = center_cex, colkey = F, ...)
  }

  invisible()
}

#' Construct the 3D tube for each curve
#'
#' This function creates a circle roughly locally perpendicular to each point along a
#' curve, where the circle has radius \code{radius}. Here, \code{len} describes how many
#' points describe this circle, where in most applications, \code{len=20} is more than enough.
#'
#' Importantly, \code{dat} needs to be in order, meaning the curve is represented by a discretization
#' of \code{nrow(dat)} points in 3-dimensional space (i.e., we require \code{ncol(dat) = 3}), whereby
#' the direct neighbors of any point (i.e., any row of \code{dat}) are the corresponding rows right before
#' and after it.
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param radius a numeric
#' @param len a numeric
#'
#' @return a list
#' @export
construct_3d_tube <- function(dat, radius, len = 20){
  stopifnot(nrow(dat) > 2, ncol(dat), is.matrix(dat))

  # remove duplicates
  dat <- .remove_duplicate_rows(dat)

  # construct all the circles
  circle_list <- .construct_all_circles(dat, radius, len = len)

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

################

#' Remove duplicated rows
#'
#' We remove rows that are both too similar to its immediate next row, as well as
#' rows that are too similar to the row after the immediate next row, both according to
#' the value in \code{tol}.
#'
#' @param dat a numeric matrix
#' @param tol a small positive number
#'
#' @return a matrix
.remove_duplicate_rows <- function(dat, tol = 1e-6){
  stopifnot(is.matrix(dat), tol > 0)

  while(TRUE){
    n <- nrow(dat)
    dist_vec1 <- sapply(1:(n-1), function(x){.l2norm(dat[x,]-dat[x+1,])})
    dist_vec2 <- sapply(1:(n-2), function(x){.l2norm(dat[x,]-dat[x+2,])})

    keep_idx <- intersect(c(1, which(dist_vec1 > tol)+1),
                          c(1, 2, which(dist_vec2 > tol)+2))

    if(length(keep_idx) == nrow(dat)) break()

    dat <- dat[keep_idx,,drop=F]
  }

  dat
}


#' Find the local direction
#'
#' Find the local direction for row \code{idx} (an index) within a dataset \code{dat}.
#' This does it by averaging the unit vector between row \code{idx} and \code{idx+1},
#' and the unit vector between row \code{idx} and \code{idx-1} (if all three indices are valid).
#' Then, it normalizes this average again, so this function always returns a unit vector
#' (i.e., a vector with L2 norm equal to 1).
#'
#' @param dat a numeric matrix
#' @param idx an index between 1 and \code{nrow(dat)}
#'
#' @return a unit vector
.find_adjacent_directions <- function(dat, idx){
  stopifnot(idx >= 1, idx <= nrow(dat))

  idx_vec <- sort(unique(pmax(pmin(idx + c(-1,0,1), nrow(dat)), 1)))
  stopifnot(length(idx_vec) %in% c(2,3))

  res <- rowMeans(sapply(1:(length(idx_vec)-1), function(i){
    vec <- dat[idx_vec[i],] - dat[idx_vec[i+1],]
    vec/.l2norm(vec)
  }))

  res/.l2norm(res)
}

#' Compute the projection
#'
#' Computes the projection of a vector \code{vec} onto the column space of \code{mat} (i.e.,
#' yields a vector with less dimensions than \code{length(vec)})
#'
#' @param vec a numeric vector
#' @param mat a numeric matrix with number of rows equal to \code{length(vec)} and with less
#' columns than rows
#'
#' @return a numeric vector of length \code{ncol(mat)}
.projection_matrix <- function(vec, mat){
  stopifnot(is.matrix(mat), length(vec) == nrow(mat), nrow(mat) >= ncol(mat))
  n <- length(vec)
  proj_mat <- diag(n) - mat %*% solve(t(mat) %*% mat) %*% t(mat)

  as.numeric(proj_mat %*% vec)
}

#' Find basis vectors that are perpendicular to a vector
#'
#' Only works in 3-dimensional space (i.e., \code{length(direction) = 3}).
#' The \code{direction} vector need not be a unit vector. However, this function
#' is guaranteed to return two unit vectors.
#'
#' @param direction a numeric vector of length 3
#'
#' @return a list of 2 unit vectors, both of length 3
.find_basis_vectors <- function(direction){
  stopifnot(length(direction) == 3)

  mat <- cbind(direction, matrix(stats::rnorm(6), ncol = 2, nrow = 3))
  for(i in 2:3) {
    mat[,i] <- .projection_matrix(mat[,i], mat[,1:(i-1), drop = F])
    mat[,i] <- mat[,i]/.l2norm(mat[,i])
  }

  list(vec1 = mat[,2], vec2 = mat[,3])
}

#' Construct a 3D circle
#'
#' Given a 3-dimensional vector \code{dat_vec} and 2 basis vectors (also 3-dimensional)
#' \code{basis_vec1} and \code{basis_vec2}, output the discretization (of \code{len} points)
#' of a circle centered at \code{dat_vec} of radius \code{radius} that lies on the 2-dimensional
#' plane formed by \code{basis_vec1} and \code{basis_vec2}.
#'
#' For this function, if \code{check = F}, we explicitly assume that \code{basis_vec1} and \code{basis_vec2}
#' are unit vectors that are perpendicular, without checking. Use this with caution
#'
#' @param dat_vec a vector of length 3
#' @param radius a positive numeric
#' @param basis_vec1 a vector of length 3
#' @param basis_vec2 a vector of length 3
#' @param len positive integer
#' @param check boolean
#'
#' @return a matrix with \code{len} rows and 3 columns
.construct_3d_circle <- function(dat_vec, radius, basis_vec1, basis_vec2,
                                 len = 20, check = T){
  stopifnot(length(dat_vec) == 3, length(basis_vec1) == 3, length(basis_vec2) == 3,
            radius > 0)

  if(check){
    basis_vec1 <- .l2norm(basis_vec1)
    basis_vec2 <- basis_vec2 - (basis_vec2 %*% basis_vec1)*basis_vec1
    basis_vec2 <- .l2norm(basis_vec2)
  }

  seq_vec <- seq(0, 2*pi, length.out = len)
  x <- radius*cos(seq_vec)
  y <- radius*sin(seq_vec)

  mat <- t(cbind(basis_vec1, basis_vec2) %*% rbind(x,y))
  for(i in 1:nrow(mat)){
    mat[i,] <- mat[i,]+dat_vec
  }

  mat
}

#' Construct all the circle for a given dataset
#'
#' The dataset \code{dat} can only have 3 columns. This is essentially
#' a helper function to streamline the usage of \code{eSVD:::.find_adjacent_directions},
#' \code{eSVD:::.find_basis_vectors} and \code{eSVD:::.construct_3d_circle}.
#'
#' @param dat a numeric matrix with 3 columns
#' @param radius a positive numeric
#' @param len a positive integer
#'
#' @return a list of matrices, each with \code{len} rows and 3 columns
.construct_all_circles <- function(dat, radius, len = 20){
  stopifnot(ncol(dat) == 3)

  n <- nrow(dat)
  lapply(1:n, function(x){
    direction <- .find_adjacent_directions(dat, x)
    res <- .find_basis_vectors(direction)
    basis_vec1 <- res$vec1; basis_vec2 <- res$vec2
    .construct_3d_circle(dat[x,], radius, basis_vec1, basis_vec2, len = len, check = F)
  })
}

#' Align two circles
#'
#' Given two circles in 3-dimensional space, represented as 2 different matrices with equal number of rows
#' and 3 columns, shuffle the row indices in \code{mat2} around so the rows of \code{mat2}
#' are (in some sense) aligned with the rows of \code{mat1}. This function has a rather
#' complicated-looking subroutine using \code{eSVD:::.subset_indices} since \code{mat1} could be
#' encoding a discretization of a circle in a clockwise direction while \code{mat2} could be
#' encoding a discretization of a circle in a counter-clockwise direction.
#'
#' The code does not explicitly check that \code{mat1} and \code{mat2} contain points that
#' are essentially an ordered discretization of a circle, but it relies on this fact.
#'
#' @param mat1 a numeric matrix with 3 columns
#' @param mat2 a numeric matrix with 3 columns and the \code{nrow(mat1)} number of rows.
#'
#' @return \code{mat2} with a new permutation of its rows
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

#' Construct a mapping from one set of indices to another
#'
#' Given an index \code{n} and two indices \code{anchor1} and \code{anchor2} that are smaller
#' than \code{n}, output a vector \code{to} and \code{from} (as a list) where
#' the first element in \code{to} is \code{anchor1} and the first element of \code{from}
#' is \code{anchor2}. Then, if \code{direction_forward=TRUE}, while the indicies of
#' \code{to} are counting up, the indicies of \code{from}
#' are also counting up. Otherwise (meaning \code{direction_forward=FALSE}),
#' while the indicies of
#' \code{to} are counting up, the indicies of \code{from}
#' are also counting down.
#'
#' In either case, all the indices here are (roughly) modulo \code{n}.
#'
#' @param n a positive integer
#' @param anchor1 a positive integer less than \code{n}
#' @param anchor2 a positive integer less than \code{n}
#' @param direction_forward a boolean
#'
#' @return a list of two vectors, each of length \code{n}, named \code{to} and \code{from}.
.subset_indices <- function(n, anchor1, anchor2, direction_forward = T){
  stopifnot(n >= anchor1, n >= anchor2, all(c(n, anchor1, anchor2) %% 1 == 0),
            all(c(n, anchor1, anchor2) > 0))

  vec1 <- anchor1:(anchor1+n-1)
  if(direction_forward){
    vec2 <- anchor2:(anchor2+n-1)
  } else {
    vec2 <- anchor2:(anchor2-n+1)
  }

  stopifnot(length(vec1) == n, length(vec2) == n)

  mat <- rbind(vec1, vec2)
  mat[mat > n] <- mat[mat > n] - n
  mat[mat < 1] <- mat[mat < 1] + n

  stopifnot(all(mat >= 1), all(mat <= n))

  list(to = mat[1,], from = mat[2,])
}







