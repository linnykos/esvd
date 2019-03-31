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
#' @param shrink shrinkage factor
#' @param thresh parameter to determine convergence
#' @param max_iter maximum number of iterations
#' @param b parameter for the kernel function (when smoothing)
#' @param upscale_vec vector of positive numbers, one for each cluster
#'
#' @return a list containing the lineages under \code{lineages},
#' the list of curves as \code{principal_curve} objects under
#' \code{curves} and the clustering matrix under \code{cluster_mat}
#' @export
slingshot <- function(dat, cluster_labels, starting_cluster,
                      cluster_group_list = NA,
                      shrink = 1, thresh = 0.001, max_iter = 15, b = 1,
                      upscale_vec = NA){
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = starting_cluster,
                            cluster_group_list = cluster_group_list)
  curves <- .get_curves(dat, cluster_labels, lineages, shrink = shrink,
                        thresh = thresh, max_iter = max_iter, b = b, upscale_vec = upscale_vec)

  list(lineages = lineages, curves = curves, cluster_mat = cluster_mat)
}

#' Estimate the slingshot curves
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param cluster_labels vector of cluster labels, where
#' the cluster labels are consecutive positive integers from 1 to
#' \code{max(cluster_labels)}
#' @param lineages output of \code{.get_lineage()}
#' @param shrink shrinkage factor
#' @param thresh parameter to determine convergence
#' @param max_iter maximum number of iterations
#' @param b parameter for the kernel function (when smoothing)
#' @param upscale_vec vector of positive numbers, one for each cluster
#'
#' @return a list of \code{principal_curve} objects
.get_curves <- function(dat, cluster_labels, lineages, shrink = 1,
                        thresh = 0.001, max_iter = 15, b = 1,
                        upscale_vec = NA){
  stopifnot(shrink >= 0 & shrink <= 1)

  if(!any(is.na(upscale_vec))){
    idx_all <- unlist(lapply(1:max(cluster_labels, na.rm = T), function(x){
      idx <- which(cluster_labels == x)
      sample(idx, upscale_vec[x]*length(idx), replace = T)
    }))
    dat <- dat[idx_all,]
    cluster_labels <- cluster_labels[idx_all]
  }

  ### setup
  num_lineage <- length(lineages)
  if(any(is.na(cluster_labels))) {
    idx <- which(is.na(cluster_labels))
    dat <- dat[-idx,]
    cluster_labels <- cluster_labels[-idx]
  }
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  cluster_vec <- 1:ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)

  W <- .initialize_weight_matrix(cluster_mat, lineages)

  ### initial curves are piecewise linear paths through the tree
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
  res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)
  pcurve_list <- res$pcurve_list; D <- res$D

  if(length(lineages) == 1) {
    s_list <- lapply(1:num_lineage, function(lin){
      sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)
      .smoother_func(pcurve_list[[lin]]$lambda, dat[sample_idx,,drop = F], b = b)
    })

    res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)
    pcurve_list <- res$pcurve_list; D <- res$D
    return(pcurve_list)
  }

  ### determine curve hierarchy
  avg_order <- .initialize_curve_hierarchy(lineages, cluster_vec)

  ### track distances between curves and data points to determine convergence
  dist_new <- sum(abs(D[W>0]))
  dist_old <- Inf

  iter <- 1
  while (abs((dist_old - dist_new) >= thresh * dist_old) && iter < max_iter){
    dist_old <- dist_new

    ### predict each dimension as a function of lambda (pseudotime)
    s_list <- lapply(1:num_lineage, function(lin){
      sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)
      .smoother_func(pcurve_list[[lin]]$lambda, dat[sample_idx,,drop = F], b = b)
    })

    res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)
    pcurve_list <- res$pcurve_list; D <- res$D
    dist_new <- sum(D[W>0], na.rm=TRUE)

    # shrink together lineages near shared clusters
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
                         pct_shrink[[i]][[j]] * shrink, dat)
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

  pcurve_list
}

##################

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
#' @param centers matrix output of \code{.compute_cluster_center()} that's
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
#' @return
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
#'
#' @return a list that contains the \code{num_lineage} curves as \code{principal_curve}
#' and a distance matrix (\code{D}) that contains the squared distance of each point
#' to its repsective lineage curve
.refine_curve_fit <- function(dat, curve_list, lineages, W, cluster_mat){
  stopifnot(nrow(dat) == nrow(W), ncol(W) == length(lineages))

  n <- nrow(dat); num_lineage <- length(lineages)
  D <- matrix(NA, nrow = n, ncol = num_lineage)
  cluster_vec <- 1:ncol(cluster_mat)

  pcurve_list <- vector("list", num_lineage)
  for(lin in 1:num_lineage){
    sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)

    pcurve <- princurve::project_to_curve(dat[sample_idx, ,drop = FALSE],
                                          s = curve_list[[lin]],
                                          stretch = 9999)
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
#' Does a kernel smoothing of \code{dat} based on the ordering given
#' by \code{lambda}. The smoothing is done via a Gaussian kernel
#' with bandwidth \code{b}.
#'
#' @param lambda vector given by one component of the \code{pricipal_curve} object
#' @param dat a \code{n} by \code{d} matrix. Here, \code{n} could be a subset of the
#' original dataset
#' @param b kernel bandwith
#'
#' @return a smoothed \code{n} by \code{d} matrix
.smoother_func <- function(lambda, dat, b = 1){
  stopifnot(length(lambda) == nrow(dat))
  ord <- order(lambda, decreasing = F)
  lambda <- lambda[ord]
  dat <- dat[ord,]

  kernel_func <- function(x, y){exp(-(x-y)^2/b)}

  t(sapply(1:length(lambda), function(i){
    weights <- kernel_func(lambda[i], lambda)
    as.numeric(colSums(diag(weights) %*% dat))/sum(weights)
  }))
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
      stats::approx(pcurve$lambda, pcurve$s[,jj], xout = lambdas_all)$y
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
  n <- nrow(pcurve$s)
  p <- ncol(pcurve$s)

  stopifnot(length(pct) == length(pcurve$lambda))
  # pct <- pct[which(pcurve$lambda_long > 0)]
  lambda <- pcurve$lambda
  s <- sapply(1:p, function(i){
    orig <- pcurve$s[,i]
    avg <- stats::approx(x = avg_curve$lambda,
                     y = avg_curve$s[,i], xout = lambda,
                     rule = 2)$y
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

#' Compute the cluster centers
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param cluster_mat a 0-1 matrix that is \code{n} by \code{k}
#'
#' @return a \code{k} by \code{d} matrix
.compute_cluster_center <- function(dat, cluster_mat){
  mat <- t(sapply(1:ncol(cluster_mat), function(x){
    idx <- which(cluster_mat[,x] == 1)
    colMeans(dat[idx,,drop=F])
  }))
  rownames(mat) <- colnames(cluster_mat)
  mat
}
