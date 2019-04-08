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

  t(cbind(basis_vec1, basis_vec2) %*% rbind(x,y))
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

  anchor1a <- which.min(res[1,])
  anchor1b <- res[2,anchor1a]

  anchor2a <- ifelse(anchor1a == ncol(res), anchor1a - 1, anchor1a + 1)
  anchor2b <- res[2,anchor2b]
  stopifnot(abs(anchor1b - anchor2b) == 1) # how to deal w/ this if false?

  # forward direction
  if(anchor2b > anchor1b){
    mat2_new <- matrix(0, ncol = ncol(mat2), nrow = nrow(mat2))
    mat2_new[anchor1a:n,] <- mat2[anchor1b:(n-anchor1a+1),]
    mat2_new[1:anchor1a,] <- mat2[anchor1a:anchor1b,]
  }
}

.subset_indices <- function(n, anchor1a, anchor1b, direction_forward = T){
  min_dist <- min(n - anchor1a + 1, n - anchor1b + 1)
  missing_bool <- (anchor1a < anchor1b)
  idx_set1 <- list(from = anchor1a:(anchor1a + min_dist - 1),
                   to = anchor1b:(anchor1b + min_dist - 1))

  if(missing_bool){
    missing_dist <- (n - anchor1a + 1) - min_dist
    stopifnot(missing_dist > 0)
    idx_set2 <- list(from = (anchor1a + min_dist):n,
                     to = 1:(missing_dist-1))
  } else {
    idx_set2 <- list(from = NA, to = NA)
  }


}





