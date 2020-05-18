#' Bootstrap the curves via resampleing of each cluster
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}
#' @param cluster_group_list list denoting the hierarchy and order of the clusters
#' @param lineages list of lineages (for example, estimated by the \code{slingshot} function)
#' @param trials numeric
#' @param ncores numeric
#' @param verbose boolean
#' @param ... additional parameters for \code{slingshot}
#'
#' @return a list
#' @export
bootstrap_curves <- function(dat, cluster_labels, cluster_group_list, lineages,
                             trials = 100, ncores = NA,
                             verbose = F, ...){
  stopifnot(!any(is.na(cluster_group_list)))

  if(!is.na(ncores)) doMC::registerDoMC(cores = ncores)

  # first do the resampling
  res <- .resample_all(dat, cluster_labels, cluster_group_list, lineages, upscale_factor)
  dat2 <- res$dat; cluster_labels <- res$cluster_labels; idx_all <- res$idx_all

  # then estimate all the different curves
  func <- function(x){
    set.seed(10*x)
    if(verbose & x %% floor(trials/10) == 0) print('*')
    dat2 <- dat
    for(i in 1:length(unique(cluster_labels))){
      idx <- which(cluster_labels == i)
      idx2 <- sample(idx, length(idx), replace = T)
      dat2[idx,] <- dat[idx2,]
    }

    .get_curves(dat2, cluster_labels, lineages = lineages, verbose = verbose, ...)
  }

  if(is.na(ncores)){
    lapply(1:trials, func)
  } else {
    x <- 0 #bookkeeping purposes
    foreach::"%dopar%"(foreach::foreach(x = 1:trials), func(x))
  }
}

#' Compute the bootstrap uncertainty
#'
#' Try to compute the "radius" of uncertainty by conditioning on the same curves
#'
#' @param target_curve_list a list (for example, the \code{curves} output from the \code{slingshot} function)
#' @param bootstrap_curve_list a list
#' (for example, the output from the \code{bootstrap_curves} function)
#' @param ncores numeric
#' @param verbose boolean
#'
#' @return a numeric
#' @export
compute_curve_sd <- function(target_curve_list, bootstrap_curve_list, ncores = NA, verbose = F){
  stopifnot(length(target_curve_list) == length(bootstrap_curve_list[[1]]),
            length(unique(sapply(bootstrap_curve_list, length))) == 1)
  num_curves <- length(target_curve_list)

  if(!is.na(ncores)) doMC::registerDoMC(cores = ncores)

  # discretize all the curves
  target_curve_list <- sapply(target_curve_list, function(curve){
    .discretize_curve_by_pseudotime(s_mat = curve$s, pseudotime_vec = curve$lambda)
  })

  for(i in 1:length(bootstrap_curve_list)){
    bootstrap_curve_list[[i]] <- sapply(bootstrap_curve_list[[i]], function(curve){
      .discretize_curve_by_pseudotime(s_mat = curve$s, pseudotime_vec = curve$lambda)
    })
  }

  # for each curve in the target curve, find the minimum distance
  mat_list <- lapply(1:num_curves, function(i){
    if(verbose) print(paste0("Starting curve ", i))
    curve_mat <- target_curve_list[[i]]$s
    curve_mat_collection <- lapply(bootstrap_curve_list, function(curve){curve$s})

    .compute_l2_curve(curve_mat, curve_mat_collection, ncores = ncores, verbose = verbose)
  })


  # output the matrix of minimum distances (position on curve x distance to bootstrap matrix)

  list(target_curve_list = target_curve_list, mat_list = mat_list)
}

#####################

.discretize_curve_by_pseudotime <- function(s_mat, pseudotime_vec, resolution = min(length(pseudotime_vec), 1000)){
  stopifnot(nrow(s_mat) == length(pseudotime_vec), nrow(s_mat) > resolution, resolution > 1)
  n <- length(pseudotime_vec)
  ord_vec <- order(pseudotime_vec, decreasing = F)
  s_mat <- s_mat[ord_vec,]
  pseudotime_vec <- pseudotime_vec[ord_vec]

  tmp_mat <- data.frame(pseudotime = pseudotime_vec, idx = 1:n)
  tmp_mat <- tmp_mat[-duplicated(tmp$pseudotime),]
  idx <- tmp_mat$idx[round(seq(1, nrow(tmp_mat), length.out = resolution))]

  list(s = s_mat[idx,], lambda = tmp_mat$pseudotime)
}

# for every point in our_mat, find its l2 distance to its closest neighbor in all curves in mat_collection
.compute_l2_curve <- function(mat, mat_collection, verbose = F, ncores = NA){
  len <- length(mat_collection)

  func <- function(mat1, mat2){
    n <- nrow(mat1)

    dist_vec <- apply(mat1, 1, function(vec1){
      min(apply(mat2, 1, function(vec2){.l2norm(vec1 - vec2)}))
    })

    dist_vec
  }

  # helper function simply for verbose purposes
  func2 <- function(mat1, i){
    if(verbose & i %% floor(len/10) == 0) cat('*')
    func(mat1, mat_collection[[i]])
  }

  if(is.na(ncores)){
    dist_mat <- matrix(NA, nrow = nrow(mat), ncol = len)
    for(i in 1:len){
      dist_mat[,i] <- func2(mat,i)
    }
  } else {
    i <- 0 #bookkeeping purposes
    dist_mat <- do.call(cbind, foreach::"%dopar%"(foreach::foreach(i = 1:len), func2(mat, i)))
  }

  dist_mat
}
