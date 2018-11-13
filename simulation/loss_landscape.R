rm(list=ls())
load("../simulation/fit.RData")

.evaluate_objective.gaussian(dat, fit$u_mat, fit$v_mat)
.evaluate_objective.gaussian(dat, res$cell_mat, res$gene_mat)

.convert_to_vector <- function(u_mat, v_mat){
  c(as.numeric(u_mat), as.numeric(v_mat))
}

.convert_to_matrix <- function(vec, n, d, k){
  stopifnot(length(vec) == n*k + d*k)
  u_mat <- matrix(vec[1:(n*k)], nrow = n, ncol = k)
  v_mat <- matrix(vec[(n*k+1):length(vec)], nrow = d, ncol = k)

  list(u_mat = u_mat, v_mat = v_mat)
}

.compute_directions <- function(list1, list2, list3,
                                vec1 = NA, vec2 = NA, vec3 = NA){
  n <- nrow(list1[[1]]); d <- nrow(list1[[2]]); k <- ncol(list1[[1]])

  if(all(is.na(vec1))) vec1 <- .convert_to_vector(list1[[1]], list1[[2]])
  if(all(is.na(vec2))) vec2 <- .convert_to_vector(list2[[1]], list2[[2]])
  if(all(is.na(vec3))) vec3 <- .convert_to_vector(list3[[1]], list3[[2]])

  # compute the directions
  dir1 <- vec2 - vec1; dist1 <- .l2norm(dir1); dir1 <- dir1/dist1

  # compute the orthogonal component
  dir2 <- vec3 - vec1
  dir2 <- as.numeric(dir2 %*% (diag(length(dir2)) - dir1 %*% t(dir1)))
  dist2 <- .l2norm(dir2); dir2 <- dir2/dist2
  stopifnot(abs(as.numeric(dir1 %*% dir2)) <= 1e-6)

  list(vec1 = vec1, dir1 = dir1, dir2 = dir2, dist1 = dist1, dist2 = dist2,
       n = n, d = d, k = k)
}

.compute_grid <- function(dat, vec1, dir1, dir2, n, d, k,
                          xrange, yrange, grid_size = 101, verbose = T){
  # now set up a coordinate system
  grid_size <- grid_size
  x_seq <- seq(xrange[1], xrange[2], length.out = grid_size)
  y_seq <- seq(yrange[1], yrange[2], length.out = grid_size)
  z_mat <- matrix(NA, grid_size, grid_size)
  rownames(z_mat) <- y_seq; colnames(z_mat) <- x_seq

  for(j in 1:grid_size){
    if(verbose & j %% floor(grid_size / 10) == 0) cat('*')

    for(i in 1:grid_size){
      tmp <- .convert_to_matrix(vec1 + x_seq[i]*dir1 + y_seq[j]*dir2, n, d, k)
      z_mat[i,j] <- tryCatch({.evaluate_objective.gaussian(dat, tmp$u_mat, tmp$v_mat)},
                             error = function(e){NA})
    }
  }

  z_mat
}

.generate_breaks <- function(z_mat, len = 100,
                             min_val = NA, max_val = NA,
                             mid_point = NA, split = 0.75){
  tmp <- z_mat
  idx <- which(!is.na(tmp))

  if(is.na(min_val)) min_val <- min(tmp[idx])
  if(is.na(max_val)) max_val <- max(tmp[idx])
  stopifnot(min_val <= min(tmp[idx]), max_val >= max(tmp[idx]), min_val < max_val)

  if(is.na(mid_point)){
    dif <- max_val - min_val
    min_val + exp(seq(log(min(dif/100, 1e-3)), log(dif), length.out = len))
  } else {
    stopifnot(min_val < mid_point, mid_point < max_val, 0 < split, split < 1)
    dif <- mid_point - min_val
    sort(c(min_val + exp(seq(log(min(dif/100, 1e-3)), log(dif), length.out = round(split*len))),
      seq(mid_point+1, max_val, length.out = round((1-split)*len))))
  }
}

##########

# try an extreme example

list1 <- list(fit$u_mat, fit$v_mat)
list2 <- list(-fit$u_mat, -fit$v_mat)
list3 <- list(res$cell_mat, res$gene_mat)
tmp <- .compute_directions(list1, list2, list3)
max_dist <- max(tmp$dist1, tmp$dist2)
z_mat <- .compute_grid(dat, vec1 = tmp$vec1, dir1 = tmp$dir1, dir2 = tmp$dir2,
                       n = tmp$n, d = tmp$d, k = tmp$k,
                       xrange = c(-0.1*max_dist, max_dist), yrange = c(-0.1*max_dist, max_dist))

image(.rotate(z_mat))


library(plot3D)
par(mar = rep(0,4))
zero <- .convert_to_matrix(rep(1e-3, tmp$n*tmp$k + tmp$d*tmp$k), tmp$n, tmp$d, tmp$k)
mid_point <- .evaluate_objective.gaussian(dat, zero$u_mat, zero$v_mat)
breaks_vec <- .generate_breaks(z_mat, len = 100, mid_point = mid_point,
                               split = 0.8)
plot3D::persp3D(z = z_mat, breaks = breaks_vec, zlim = range(breaks_vec))

##############

list1 <- list(fit$u_mat, fit$v_mat)
list2 <- list(fit$res_list[[9]]$u_mat, fit$res_list[[9]]$v_mat)
list3 <- list(fit$res_list[[8]]$u_mat, fit$res_list[[8]]$v_mat)
tmp <- .compute_directions(list1, list2, list3)
max_dist <- max(tmp$dist1, tmp$dist2)
z_mat <- .compute_grid(dat, vec1 = tmp$vec1, dir1 = tmp$dir1, dir2 = tmp$dir2,
                       n = tmp$n, d = tmp$d, k = tmp$k,
                       xrange = c(-0.1*max_dist, max_dist), yrange = c(-0.1*max_dist, max_dist))
plot3D::persp3D(z = z_mat, breaks = breaks_vec)

###################

# try a random direction
list1 <- list(fit$u_mat, fit$v_mat)
vec1 <- .convert_to_vector(list1[[1]], list1[[2]])
n <- nrow(list1[[1]]); d <- nrow(list1[[2]]); k <- ncol(list1[[1]])

iter <- 1
min_perc <- Inf

while(TRUE){
  if(iter %% 20 == 0) print(paste0("Iteration: ", iter, ", Percentage: ", min_perc))
  set.seed(10*iter)
  dir1 <- rnorm(n*k+d*k); dir1 <- dir1/.l2norm(dir1)
  dir2 <- rnorm(n*k+d*k); dir2 <- dir2/.l2norm(dir2)
  dir2 <- as.numeric(dir2 %*% (diag(length(dir2)) - dir1 %*% t(dir1)))
  dist2 <- .l2norm(dir2); dir2 <- dir2/dist2
  max_dist <- 5
  z_mat <- .compute_grid(dat, vec1 = vec1, dir1 = dir1, dir2 = dir2,
                         n = tmp$n, d = tmp$d, k = tmp$k,
                         xrange = c(-max_dist, max_dist),
                         yrange = c(-max_dist, max_dist), grid_size = 5,
                         verbose = F)
  perc <- length(which(is.na(z_mat)))/prod(dim(z_mat))
  if(perc <= min_perc) min_perc <- perc
  if(perc <= 0.1) break()
  iter <- iter + 1
}


