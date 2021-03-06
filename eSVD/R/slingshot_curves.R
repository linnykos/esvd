# code adapted from https://github.com/kstreet13/slingshot

#' Use slingshot to estimate the cell development trajectories
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}
#' @param starting_cluster the "origin" cluster that all the lineages will start
#' from
#' @param cluster_group_list list denoting the hierarchy and order of the clusters
#' @param squared boolean on whether or not to square the distance matrix
#' @param shrink shrinkage factor
#' @param stretch strecth factor
#' @param thresh parameter to determine convergence
#' @param max_iter maximum number of iterations
#' @param upscale_factor positive numeric (between 0 and 1) that controls how much to upweight the clusters,
#' with 1 being (almost) equal cluster sizes and 0 being no upweighting. This does not affect the
#' estimation of the lineages via \code{.get_lineages}, only the estimation of the curves
#' given the lineages via \code{.get_curves}.
#' @param verbose boolean
#'
#' @return a list containing the lineages under \code{lineages},
#' the list of curves as \code{principal_curve} objects under
#' \code{curves} and the clustering matrix under \code{cluster_mat}
#' @export
slingshot <- function(dat, cluster_labels, starting_cluster,
                      cluster_group_list = NA,
                      squared = F, shrink = 1, stretch = 9999,
                      thresh = 0.001, max_iter = 15,
                      upscale_factor = NA, verbose = F){
  stopifnot(ncol(dat) >= 2)

  if(verbose) print("Starting to infer lineages")
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = starting_cluster,
                            cluster_group_list = cluster_group_list,
                            squared = squared)

  # resample matrix
  res <- .resample_all(dat, cluster_labels, cluster_group_list, lineages, upscale_factor)
  dat2 <- res$dat; cluster_labels <- res$cluster_labels; idx_all <- res$idx_all

  if(verbose) print("Starting to infer curves")
  res <- .get_curves(dat2, cluster_labels, lineages,
                     shrink = shrink, stretch = stretch, thresh = thresh, max_iter = max_iter,
                     verbose = verbose)
  curves <- res$pcurve_list

  structure(list(lineages = lineages, curves = curves, idx = idx_all), class = "slingshot")
}

#' Prepare trajectories for usage
#'
#' This function only extracts the \code{s} component of the curve, where the curves
#' are \code{principal_curve} objects, typically returned within the output of \code{eSVD::slingshot}.
#' This function truncates the curve (starting from the pseudotime index of 0, according
#' to the \code{ord} component of \code{curve}) so the total length is less that \code{target_length}
#'
#' @param curve \code{principal_curve} object
#' @param target_length target length of the curve.
#'
#' @return numeric matrix
#' @export
prepare_trajectory <- function(curve, target_length){
  stopifnot(class(curve) == "principal_curve")

  s_mat <- curve$s[curve$ord,]
  len_cumulative <- cumsum(sapply(2:nrow(s_mat), function(i){.l2norm(s_mat[i-1,] - s_mat[i,])}))

  if(max(len_cumulative) <= target_length){
    return(s_mat)
  }

  keep_idx <- max(which(len_cumulative <= target_length))
  s_mat[1:keep_idx,]
}

#' Compute the cluster centers
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param cluster_mat a 0-1 matrix that is \code{n} by \code{k}
#'
#' @return a \code{k} by \code{d} matrix
#' @export
compute_cluster_center <- function(dat, cluster_mat){
  mat <- t(sapply(1:ncol(cluster_mat), function(x){
    idx <- which(cluster_mat[,x] == 1)
    colMeans(dat[idx,,drop=F])
  }))
  rownames(mat) <- colnames(cluster_mat)
  mat
}



###################################################

#' Estimate the slingshot curves
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}
#' @param lineages output of \code{.get_lineage()}
#' @param shrink shrinkage factor
#' @param stretch stretch factor
#' @param thresh parameter to determine convergence
#' @param max_iter maximum number of iterations
#' @param verbose boolean
#'
#' @return a list of \code{principal_curve} objects
.get_curves <- function(dat, cluster_labels, lineages,
                        shrink = 1, stretch = 9999,
                        thresh = 0.001, max_iter = 15,
                        verbose = F){
  ### setup
  num_lineage <- length(lineages)
  names(lineages) <- paste0("Lineage", 1:num_lineage)
  if(any(is.na(cluster_labels))) {
    idx <- which(is.na(cluster_labels))
    dat <- dat[-idx,]
    cluster_labels <- cluster_labels[-idx]
  }
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  cluster_vec <- 1:ncol(cluster_mat)
  centers <- compute_cluster_center(dat, cluster_mat)

  W <- .initialize_weight_matrix(cluster_mat, lineages)

  ### initial curves are piecewise linear paths through the tree
  if(verbose) print("Starting to initialize curves")
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
  res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat, stretch)
  pcurve_list <- res$pcurve_list; D <- res$D

  if(length(lineages) == 1) {
    s_list <- lapply(1:num_lineage, function(lin){
      sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)
      .smoother_func(pcurve_list[[lin]]$lambda, dat[sample_idx,,drop = F], verbose = verbose)
    })

    res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat, stretch)
    pcurve_list <- res$pcurve_list; D <- res$D

    names(pcurve_list) <- paste('Curve', 1:length(pcurve_list), sep='')

    return(list(pcurve_list = pcurve_list))
  }

  ### determine curve hierarchy
  avg_order <- .initialize_curve_hierarchy(lineages, cluster_vec)

  ### track distances between curves and data points to determine convergence
  dist_new <- sum(abs(D[W>0]))
  dist_old <- Inf

  iter <- 1
  while (abs((dist_old - dist_new) >= thresh * dist_old) && iter < max_iter){
    if(verbose) print(paste0("On iteration ", iter))
    dist_old <- dist_new

    ### predict each dimension as a function of lambda (pseudotime)
    s_list <- lapply(1:num_lineage, function(lin){
      if(verbose) print(paste0("Smoothing lineage ", lin))
      sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)
      .smoother_func(pcurve_list[[lin]]$lambda, dat[sample_idx,,drop = F], verbose = verbose)
    })

    if(verbose) print("Refining curves")
    res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat, stretch)
    pcurve_list <- res$pcurve_list; D <- res$D
    dist_new <- sum(D[W>0], na.rm=TRUE)

    # shrink together lineages near shared clusters
    if(verbose) print("Shrinking curves together")
    if(shrink > 0){
      avg_curve_list <- vector("list", length(avg_order))
      names(avg_curve_list) <- paste0("Average", 1:length(avg_order))
      pct_shrink <- vector("list", length(avg_order))

      ### determine average curves and amount of shrinkage
      for (i in 1:length(avg_order)) {
        to_avg_curves <- lapply(avg_order[[i]], function(n_element){
          if(grepl('Lineage', n_element)){
            pcurve_list[[n_element]]
          } else {
            avg_curve_list[[n_element]]
          }
        })

        avg <- .construct_average_curve(to_avg_curves, dat)
        avg_curve_list[[i]] <- avg

        # find the indicies shared by all the curves
        common_ind <- which(rowMeans(sapply(to_avg_curves, function(crv){ crv$W > 0 })) == 1)
        pct_shrink[[i]] <- lapply(to_avg_curves, function(curve) {
          .percent_shrinkage(curve, common_ind)
        })

        pct_shrink[[i]] <- .check_shrinkage(pct_shrink[[i]])
      }

      ### do the shrinking in reverse order
      for(i in rev(1:length(avg_curve_list))){
        avg_curve <- avg_curve_list[[i]]
        to_shrink_list <-  lapply(avg_order[[i]], function(n_element){
          if(grepl('Lineage', n_element)){
            pcurve_list[[n_element]]
          } else {
            avg_curve_list[[n_element]]
          }
        })

        shrunk_list <- lapply(1:length(to_shrink_list),function(j){
          pcurve <- to_shrink_list[[j]]
          .shrink_to_avg(pcurve, avg_curve,
                         pmin(pct_shrink[[i]][[j]] * shrink, 1), dat)
        })

        for(j in 1:length(avg_order[[i]])){
          ns <- avg_order[[i]][[j]]
          if(grepl('Lineage', ns)){
            pcurve_list[[ns]] <- shrunk_list[[j]]
          } else {
            avg_curve_list[[ns]] <- shrunk_list[[j]]
          }
        }
      }
    }

    iter <- iter + 1
  }

  names(pcurve_list) <- paste('Curve', 1:length(pcurve_list), sep='')

  list(pcurve_list = pcurve_list)
}

##################

#' Intersect lineages with cluster group list
#'
#' Recall: \code{lineages} is a list of vectors, where each vector represents an ordered
#' set of indices. \code{cluster_group_list} is an ordered list of vectors, where each vector represent
#' indices "of the same tier." This function intersects the two by finding a list of vectors where
#' each of the output set of vectors is a self-contained vector in both \code{lineages} and \code{cluster_group_list}.
#'
#' @param lineages a list of vectors
#' @param cluster_group_list a list of vectors
#'
#' @return a list of vectors
.intersect_lineages_cluster_group_list <- function(lineages, cluster_group_list){
  len <- length(cluster_group_list)
  lineage_len <- length(lineages)

  tmp <- lapply(cluster_group_list, function(cluster_group){
    if(length(cluster_group) == 1) return(list(cluster_group))

    # for each element in x, determine which lineages it is in
    lineage_char <- sapply(cluster_group, function(clust){
      idx <- which(sapply(lineages, function(lineage){clust %in% lineage}))
      paste0(idx, collapse = "-")
    })
    lineage_char <- factor(lineage_char)

    lis <- lapply(levels(lineage_char), function(lvl){
      cluster_group[which(lineage_char == lvl)]
    })

    lis
  })

  # reformat list so it's just one big list, not nested lists
  .flatten_list(tmp)
}

#' Flatten a list of lists
#'
#' Code from \url{https://stackoverflow.com/questions/8139677/how-to-flatten-a-list-to-a-list-without-coercion/8139959#8139959}
#'
#' @param x a list of lists
#'
#' @return a list
.flatten_list <- function(x) {
  stopifnot(is.list(x))
  stopifnot(all(sapply(x, class) == "list"))

  len <- sum(rapply(x, function(x) 1L))
  y <- vector('list', len)
  i <- 0L
  rapply(x, function(x) { i <<- i+1L; y[[i]] <<- x })
  y
}

#' Resampling rows in a dattaset according to the amount of upscaling
#'
#' This is the primary function used in \code{eSVD::slingshot} to decide which
#' rows in \code{dat} to resample (primarily useful for getting better-looking curves
#' in \code{eSVD:::.get_curves})
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}
#' @param cluster_group_list list denoting the hierarchy and order of the clusters
#' @param lineages list of vectors containing cluster labels
#' @param upscale_factor numeric scalar
#'
#' @return a list containing a new version of \code{dat} and \code{cluster_labels} but also
#' \code{idx_all} which stores the row-indices of the original \code{dat} used
.resample_all <- function(dat, cluster_labels, cluster_group_list, lineages, upscale_factor){
  if(all(is.na(cluster_group_list))) cluster_group_list <- list(sort(unique(cluster_labels)))

  if(!any(is.na(upscale_factor))){
    # intersect lineages with cluster_group_list to determine group sizes
    cluster_intersection <- .intersect_lineages_cluster_group_list(lineages, cluster_group_list)
    upscale_vec <- .compute_upscale_factor(cluster_labels, cluster_intersection, upscale_factor)

    idx_all <- .construct_resample_idx(cluster_labels, upscale_vec)
    dat <- dat[idx_all,]
    cluster_labels <- cluster_labels[idx_all]
  } else {
    idx_all <- 1:nrow(dat)
  }

  list(dat = dat, cluster_labels = cluster_labels, idx_all = idx_all)
}

#' Construct resampled cluster labels
#'
#' This function is designed so all the indices in 1 through \code{length(cluster_labels)} is
#' represented at least once. Typically, all the values in \code{upscale_vec} are 1 or larger.
#'
#' @param cluster_labels  vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}
#' @param upscale_vec vector of positive values with length equal to \code{max(cluster_labels)}
#'
#' @return a vector of cluster labels (i.e. indices)
.construct_resample_idx <- function(cluster_labels, upscale_vec){
  stopifnot(max(cluster_labels) == length(upscale_vec))

  unlist(lapply(1:max(cluster_labels), function(x){
    idx <- which(cluster_labels == x)

    if(upscale_vec[x] >= 1){
      idx <-  c(idx, sample(idx, max(c(round((upscale_vec[x]-1)*length(idx)), 0)), replace = T))
    }

    idx
  }))
}

#' Compute the upscale factor
#'
#' The upscale factor is the fraction of "maximum size" over "size of cluster", raised to \code{upscale_factor}.
#' Hence, if \code{upscale_factor=1}, all the resulting clusters in \code{cluster_intersection} are of the same size.
#'
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}
#' @param cluster_intersection list of vectors, where each vector contains entries in \code{cluster_labels}
#' @param upscale_factor numeric
#'
#' @return a numeric vector equal in length to \code{length(unique(cluster_labels))}
.compute_upscale_factor <- function(cluster_labels, cluster_intersection, upscale_factor = 0.5){
  stopifnot(all(unlist(cluster_intersection) %in% cluster_labels))
  stopifnot(length(unlist(cluster_intersection)) == max(cluster_labels))

  size_vec <- sapply(cluster_intersection, function(x){length(which(cluster_labels %in% x))})
  upscale_vec <- rep(NA, length(unique(cluster_labels)))

  for(i in 1:length(cluster_intersection)){
    upscale_vec[cluster_intersection[[i]]] <- (max(size_vec)/size_vec[i])^(upscale_factor)
  }

  upscale_vec
}


#' Initialize the weight matrix
#'
#' @param cluster_mat a 0-1 matrix with n rows and k columns
#' @param lineages output of \code{.get_lineages()}
#'
#' @return 0-1 matrix of size n by \code{length(lineages)}
.initialize_weight_matrix <- function(cluster_mat, lineages){
  num_lineage <- length(lineages)
  W <- sapply(1:num_lineage, function(i){
    rowSums(cluster_mat[, lineages[[i]], drop = FALSE])
  })
  colnames(W) <- names(lineages)

  W
}

#' Constructs the order to average the lineages
#'
#' As stated in the paper, "Average curves are constructed  in a recursive
#' manner, from the latest branching events to the earliest."
#'
#' @param lineages output of \code{.get_lineages()}
#' @param cluster_vec vector of cluster labels
#'
#' @return list of lineages to average
.initialize_curve_hierarchy <- function(lineages, cluster_vec){
  C <- as.matrix(sapply(lineages, function(lin) {
    sapply(cluster_vec, function(i) {
      as.numeric(i %in% lin)
    })}))
  rownames(C) <- cluster_vec

  segmnts <- unique(C[rowSums(C)>1,,drop = FALSE])
  segmnts <- segmnts[order(rowSums(segmnts),decreasing = FALSE), ,
                     drop = FALSE]
  avg_order <- list()

  for(i in 1:nrow(segmnts)){
    idx <- which(segmnts[i,] == 1)
    avg_order[[i]] <- colnames(segmnts)[idx]
    new_col <- rowMeans(segmnts[,idx, drop = FALSE])
    segmnts <- cbind(segmnts[, -idx, drop = FALSE], new_col)
    colnames(segmnts)[ncol(segmnts)] <- paste('Average',i,sep='')
  }

  avg_order
}


#' Initial curve fit
#'
#' @param lineages list output of \code{.get_lineages()}
#' @param cluster_vec a vector from 1 to \code{k}
#' @param centers matrix output of \code{compute_cluster_center()} that's
#' \code{k} by \code{d}
#'
#' @return a list of points
.initial_curve_fit <- function(lineages, cluster_vec, centers){
  lapply(lineages, function(lineage){
    line_initial <- centers[cluster_vec %in% lineage, , drop = FALSE]
    line_initial[match(lineage, rownames(line_initial)),,drop = FALSE]
  })
}


#' Determine which samples are part of a lineage
#'
#' @param lineage one specific element of the list output of \code{.get_lineages()}
#' @param cluster_mat 0-1 matrix output of \code{.construct_cluster_matrix()} that
#' has \code{n} rows and \code{k} column
#'
#' @return a vector of numbers
.determine_idx_lineage <- function(lineage, cluster_mat){
  sort(unique(unlist(lapply(as.numeric(lineage), function(x){
    which(cluster_mat[,x] == 1)
  }))))
}

.clean_curve <- function(pcurve, W_vec, sample_idx = NA){
  stopifnot(class(pcurve) == "principal_curve")
  stopifnot(!is.matrix(W_vec))
  stopifnot(length(which(W_vec > 0)) == length(pcurve$lambda))

  # force non-negative distances
  pcurve$dist_ind <- abs(pcurve$dist_ind)
  # force pseudotime to start at 0
  pcurve$lambda <- pcurve$lambda - min(pcurve$lambda, na.rm=TRUE)
  pcurve$lambda_long <- rep(0, length(W_vec))
  pcurve$lambda_long[which(W_vec > 0)] <- pcurve$lambda
  pcurve$W <- W_vec
  pcurve$idx <- sample_idx

  pcurve
}


#' Initial fit of the curves
#'
#' This method is heavily dependent on \code{princurve::project_to_curve()}.
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param curve_list sequence of points, input for \code{princurve::project_to_curve()}
#' @param lineages list output of \code{.get_lineages()}
#' @param W a weight matrix that is \code{n} by \code{num_lineage} (number of lineages)
#' @param cluster_mat 0-1 matrix output of \code{.construct_cluster_matrix()} that
#' has \code{n} rows and \code{k} column
#' @param stretch stretch factor
#'
#' @return a list that contains the \code{num_lineage} curves as \code{principal_curve}
#' and a distance matrix (\code{D}) that contains the squared distance of each point
#' to its repsective lineage curve
.refine_curve_fit <- function(dat, curve_list, lineages, W, cluster_mat, stretch = 9999){
  stopifnot(nrow(dat) == nrow(W), ncol(W) == length(lineages))

  n <- nrow(dat); num_lineage <- length(lineages)
  D <- matrix(NA, nrow = n, ncol = num_lineage)
  cluster_vec <- 1:ncol(cluster_mat)

  pcurve_list <- vector("list", num_lineage)
  for(lin in 1:num_lineage){
    sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)

    pcurve <- princurve::project_to_curve(dat[sample_idx, ,drop = FALSE],
                                          s = curve_list[[lin]],
                                          stretch = stretch)
      # note: princurve::project_to_curve changes the input s
    pcurve <- .clean_curve(pcurve, W[, lin], sample_idx)
    pcurve_list[[lin]] <- pcurve

    D[sample_idx,lin] <- abs(pcurve$dist_ind)
  }

  names(pcurve_list) <- paste0("Lineage", 1:num_lineage)

  list(pcurve_list = pcurve_list, D = D)
}


#' Smooth approximation based on lambdas
#'
#' Does nonparametric regression (local polynomial regression)
#' of \code{dat} column-by-column based on the ordering given
#' by \code{lambda}.
#'
#' This function re-orders the rows of \code{dat} according to \code{lambda}. While it's
#' not necessary for this function, other functions that depend on this function's output
#' will require this.
#'
#' @param lambda vector given by one component of the \code{principal_curve} object
#' @param dat a \code{n} by \code{d} matrix. Here, \code{n} could be a subset of the
#' original dataset
#' @param verbose boolean
#'
#' @return a smoothed \code{n} by \code{d} matrix that is ordered by \code{lambda}
.smoother_func <- function(lambda, dat, verbose = F){
  stopifnot(length(lambda) == nrow(dat))

  # order the data
  ord <- order(lambda, decreasing = F)
  lambda <- lambda[ord]
  dat <- dat[ord,]

  sapply(1:ncol(dat), function(j){
    tmp_df <- data.frame(y = dat[,j], x = lambda)
    res <- stats::loess(y ~ x, data = tmp_df, span = 0.5)
    res$fitted
  })
}

#' Construct average curve
#'
#' @param pcurve_list list of outputs of \code{princurve::project_to_curve()}
#' @param dat a \code{n} by \code{d} matrix
#'
#' @return a \code{principal_curve} object
.construct_average_curve <- function(pcurve_list, dat){
  n <- nrow(pcurve_list[[1]]$s)
  p <- ncol(pcurve_list[[1]]$s)
  lambdas_all <- unique(unlist(lapply(pcurve_list, function(pcv){pcv$lambda})))
  max_shared_lambda <- min(sapply(pcurve_list, function(pcv){max(pcv$lambda)}))
  lambdas_all <- sort(lambdas_all[lambdas_all <= max_shared_lambda])

  # interpolate all the curves so they're parameterized on the same points
  pcurves_dense <- lapply(pcurve_list, function(pcurve){
    sapply(1:p, function(jj){
      suppressWarnings(stats::approx(pcurve$lambda, pcurve$s[,jj], xout = lambdas_all)$y)
    })
  })

  avg <- sapply(1:p, function(j){
    dim_all <- sapply(1:length(pcurves_dense),function(i){
      pcurves_dense[[i]][,j]
    })
    rowMeans(dim_all)
  })

  W_vec <- rowMeans(sapply(pcurve_list, function(p){ p$W }))
  sample_idx <- which(W_vec > 0)
  avg_curve <- princurve::project_to_curve(dat[sample_idx,,drop=F], avg)
  avg_curve <- .clean_curve(avg_curve, W_vec, sample_idx)

  avg_curve
}

#' Determine the percentage shrinkage (a non-decreasing function)
#'
#' Determine the lambda's to use for interpolation based on the IQR for
#' lambda at \code{common_idx} indicies. Then use a linear interpolation
#' via \code{approx} to assign all the lambdas (even those not in \code{common_idx})
#' a relevant value based on the CDF of the cosine kernel.
#'
#' @param pcurve output of \code{princurve::project_to_curve()}
#' @param common_idx indices to use
#'
#' @return vector
.percent_shrinkage <- function(pcurve, common_idx){
  lambda <- pcurve$lambda_long

  dens <- stats::density(0, bw = 1, kernel = "cosine")
  surv <- list(x = dens$x, y = (sum(dens$y) - cumsum(dens$y))/sum(dens$y))
  box_vals <- graphics::boxplot(lambda[common_idx], plot = FALSE)$stats
  if(box_vals[1] == box_vals[5]){
    pct_l <- rep(0, length(pcurve$lambda))
  } else {
    surv$x <- .scale_vector(surv$x, lower = box_vals[1], upper = box_vals[5])
    pct_l <- stats::approx(surv$x, surv$y, pcurve$lambda, rule = 2)$y
  }

  pct_l
}

#' Scale a vector to lie between two endpoints
#'
#' Takes \code{x} and scales it so the relative distances between points
#' are perseved and it ranges between \code{lower} and \code{upper}
#'
#' @param x vector
#' @param lower numeric
#' @param upper numeric
#'
#' @return vector
.scale_vector <- function(x, lower=0, upper=1){
  stopifnot(lower < upper)
  ((x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x, na.rm=TRUE)))*(upper-lower)+lower
}

.check_shrinkage <- function(pct_shrink_list){
  all_zero <- sapply(pct_shrink_list, function(pij){
    all(pij == 0)
  })
  if(any(all_zero)){
    pct_shrink_list <- lapply(pct_shrink_list, function(pij){
                                pij[] <- 0; pij
                              })
  }

  pct_shrink_list
}

#' Shrink curve to average
#'
#' @param pcurve output of \code{princurve::project_to_curve()}
#' @param avg_curve the average curve constructed by one of the elements
#' in the output of \code{.construct_average_curve()}
#' @param pct one of the elements in the output of \code{.percent_shrinkage()}
#' @param dat a \code{n} by \code{d} matrix
#'
#' @return a \code{principal_curve} object
.shrink_to_avg <- function(pcurve, avg_curve, pct, dat){
  stopifnot(all(pct >= 0), all(pct <= 1))
  n <- nrow(pcurve$s)
  p <- ncol(pcurve$s)

  stopifnot(length(pct) == length(pcurve$lambda))
  # pct <- pct[which(pcurve$lambda_long > 0)]
  lambda <- pcurve$lambda

  s <- sapply(1:p, function(i){
    orig <- pcurve$s[,i]

    # suppressing warnings to avoid errors to non-unique values to avg_curve$lambda
    suppressWarnings(avg <- stats::approx(x = avg_curve$lambda,
                     y = avg_curve$s[,i], xout = lambda,
                     rule = 2)$y)
    avg * pct + orig * (1-pct)
  })

  W_vec <- pcurve$W
  sample_idx <- pcurve$idx
  pcurve <- princurve::project_to_curve(dat[sample_idx,,drop = F],
                                        as.matrix(s[pcurve$ord,,drop = FALSE]))
  pcurve <- .clean_curve(pcurve, W_vec, sample_idx)

  pcurve
}

#' Construst cluster matrix from cluster labels
#'
#' @param cluster_labels  vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels, na.rm = T)}. Can include \code{NA}
#'
#' @return A 0-1 matrix with \code{length(cluster_labels)} rows
#' and \code{max(cluster_labels)} columns
.construct_cluster_matrix <- function(cluster_labels){
  idx <- !is.na(cluster_labels)
  stopifnot(all(cluster_labels[idx] > 0), all(cluster_labels[idx] %% 1 == 0))
  stopifnot(length(unique(cluster_labels[idx])) == max(cluster_labels[idx]))

  k <- max(cluster_labels, na.rm = T)
  n <- length(cluster_labels)

  mat <- sapply(1:k, function(x){
    tmp <- rep(0, n)
    tmp[which(cluster_labels == x)] <- 1
    tmp
  })

  colnames(mat) <- seq_len(k)
  mat
}

