bootstrap_curves <- function(dat, cluster_labels, starting_cluster,
                             cluster_group_list = NA, trials = 100, ...){
  func <- function(x){
    reduction_factor <- max(apply(dat, 2, function(x){diff(range(x))}))*.25
    set.seed(10*x)
    dat2 <- dat
    for(i in 1:length(unique(cluster_labels))){
      idx <- which(cluster_labels == i)
      idx2 <- sample(idx, length(idx), replace = T)
      dat2[idx,] <- dat[idx2,]
    }

    slingshot(dat2/reduction_factor, cluster_labels, starting_cluster = starting_cluster,
              cluster_group_list = cluster_group_list, ...)
  }

  lapply(1:trials, func)
}

compute_curve_sd <- function(target_curve_list, bootstrap_curve_list){
  num_curves <- length(target_curve_list$lineages)

  mat_list <- lapply(1:num_curves, function(i){
    print(paste0("Starting curve ", i))
    curve_mat <- target_curve_list$curves[[i]]$s[target_curve_list$curves[[i]]$ord,]
    curve_mat_collection <- .capture_curves(paste0(target_curve_list$lineages[[i]], collapse = "-"), bootstrap_curve_list)

    .compute_l2_curve(curve_mat, curve_mat_collection)
  })

  sd_vec <- sapply(mat_list, function(x){median(apply(x, 2, max))})

  list(sd_val = max(sd_vec), mat_list = mat_list)
}

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

#####################

# try to compute the "radius" of uncertainty by conditioning on the same curves
# first store all the relevant curves
.capture_curves <- function(string, curve_list){
  res <- lapply(1:length(curve_list), function(i){
    string_vec <- sapply(curve_list[[i]]$lineages, function(x){paste0(x, collapse="-")})
    idx <- which(string_vec == string)
    if(length(idx) == 0) return(NA)
    curve_list[[i]]$curves[[idx]]$s[curve_list[[i]]$curves[[idx]]$ord,]
  })

  res[which(sapply(res, length) > 1)]
}

# for every point in our_mat, find its l2 distance to its closest neighbor in all curves in our_mat_collection
.compute_l2_curve <- function(mat, mat_collection){
  n <- nrow(mat); k <- length(mat_collection)
  sapply(1:n, function(x){
    if(x %% floor(n/10) == 0) cat('*')

    vec <- mat[x,]
    sapply(1:k, function(y){
      dist_vec <- apply(mat_collection[[k]], 1, function(z){
        .l2norm(z-vec)
      })
      min(dist_vec)
    })
  })
}

#####


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





