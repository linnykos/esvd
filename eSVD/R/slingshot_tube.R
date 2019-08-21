#' Bootstrap the curves via resampleing of each cluster
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}
#' @param starting_cluster the "origin" cluster that all the lineages will start
#' from
#' @param cluster_group_list list denoting the hierarchy and order of the clusters
#' @param trials numeric
#' @param cores numeric
#' @param ... additional parameters for \code{slingshot}
#'
#' @return a list
#' @export
bootstrap_curves <- function(dat, cluster_labels, starting_cluster,
                             cluster_group_list, trials = 100, cores = NA,
                             ...){
  stopifnot(!any(is.na(cluster_group_list)))

  if(!is.na(cores)) doMC::registerDoMC(cores = cores)

  func <- function(x){
    set.seed(10*x)
    dat2 <- dat
    for(i in 1:length(unique(cluster_labels))){
      idx <- which(cluster_labels == i)
      idx2 <- sample(idx, length(idx), replace = T)
      dat2[idx,] <- dat[idx2,]
    }

    slingshot(dat2, cluster_labels, starting_cluster = starting_cluster,
              cluster_group_list = cluster_group_list, ...)
  }

  if(is.na(cores)){
    lapply(1:trials, func)
  } else {
    x <- 0 #bookkeeping purposes
    foreach::"%dopar%"(foreach::foreach(x = 1:trials), func(x))
  }
}

#' Compute the bootstrap uncertainty
#'
#' @param target_curve_list a list
#' @param bootstrap_curve_list a list
#' @param cores numeric
#'
#' @return a numeric
#' @export
compute_curve_sd <- function(target_curve_list, bootstrap_curve_list, cores = NA){
  num_curves <- length(target_curve_list$lineages)

  if(!is.na(cores)) doMC::registerDoMC(cores = cores)

  func <- function(i){
    print(paste0("Starting curve ", i))
    curve_mat <- target_curve_list$curves[[i]]$s[target_curve_list$curves[[i]]$ord,]
    curve_mat_collection <- .capture_curves(paste0(target_curve_list$lineages[[i]], collapse = "-"), bootstrap_curve_list)

    .compute_l2_curve(curve_mat, curve_mat_collection)
  }

  if(is.na(cores)){
    mat_list <- lapply(1:num_curves, func)
  } else {
    i <- 0 #bookkeeping purposes
    mat_list <- foreach::"%dopar%"(foreach::foreach(i = 1:num_curves), func(i))
  }

  sd_vec <- sapply(mat_list, function(x){
    stats::quantile(apply(x, 2, function(x){
      stats::quantile(x,probs = 0.95)
    }), probs = 0.95)
  })

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
.compute_l2_curve <- function(mat, mat_collection, verbose = F){
  n <- nrow(mat); k <- length(mat_collection)
  sapply(1:n, function(x){
    if(verbose & x %% floor(n/10) == 0) cat('*')

    vec <- mat[x,]
    sapply(1:k, function(y){
      dist_vec <- apply(mat_collection[[k]], 1, function(z){
        .l2norm(z-vec)
      })

      min(dist_vec)
    })
  })
}
