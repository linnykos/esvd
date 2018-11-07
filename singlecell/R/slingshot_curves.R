# code adapted from https://github.com/kstreet13/slingshot

#' Title
#'
#' @param dat
#' @param cluster_labels
#' @param lineages
#' @param shrink
#' @param extend
#' @param thresh parameter to determine convergence
#' @param max_iter
#' @param allow.breaks
#' @param b parameter for the kernel function (when smoothing)
#'
#' @return
#' @export
#'
#' @examples
.get_curves <- function(dat, cluster_labels, lineages, shrink = 1,
                        extend = 'y',
                        thresh = 0.001, max_iter = 15,
                        allow.breaks = TRUE, b = 1){
  stopifnot(shrink >= 0 & shrink <= 1)

  ### setup
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  cluster_vec <- 1:ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)

  W <- .initialize_weight_matrix(cluster_mat, lineages)

  ### determine curve hierarchy
  avg_order <- .initialize_curve_hierarchy(lineages, cluster_vec)

  ### initial curves are piecewise linear paths through the tree
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
  res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)
  pcurve_list <- res$pcurve_list; D <- res$D

  ### track distances between curves and data points to determine convergence
  dist_new <- sum(abs(D[W>0]))
  dist_old <- Inf

  iter <- 1
  while (abs((dist_old - dist_new) >= thresh * dist_old) && iter < max_iter){
    dist_old <- dist_new

    ### predict each dimension as a function of lambda (pseudotime)
    s_list <- lapply(1:num_lineage, function(lin){
      sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)
      .smoother_func(pcurve_list[[lin]]$lambda, dat[sample_idx,], b = b)
    })

    res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)
    pcurve_list <- res$pcurve_list; D <- res$D

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
        common_ind <- which(rowMeans(sapply(to_avg_curves, function(crv){ crv$w > 0 })) == 1)
        pct_shrink[[i]] <- lapply(to_avg, function(curve) {
          .percent_shrinkage(curve, common_ind)
        })


      }

      ### do the shrinking in reverse order
      for(j in rev(seq_along(avg_lines))){
        ns <- avg_order[[j]]
        avg <- avg_lines[[j]]
        to_shrink <- lapply(ns, function(n_element){
          if(grepl('Lineage', n_element)){
            l_ind <- as.numeric(gsub('Lineage','',n_element))
            return(pcurves[[l_ind]])
          }
          if(grepl('average',n)){
            a_ind <- as.numeric(gsub('average','',n_element))
            return(avg_lines[[a_ind]])
          }
        })

        shrunk <- lapply(seq_along(ns),function(jj){
          crv <- to_shrink[[jj]]
          .shrink_to_avg(crv, avg,
                         pct_shrink[[j]][[jj]] * shrink,
                         dat, stretch = stretch)
        })

        for(jj in seq_along(ns)){
          n <- ns[jj]
          if(grepl('Lineage',n)){
            l_ind <- as.numeric(gsub('Lineage','',n))
            pcurves[[l_ind]] <- shrunk[[jj]]
          }
          if(grepl('average',n)){
            a_ind <- as.numeric(gsub('average','',n))
            avg_lines[[a_ind]] <- shrunk[[jj]]
          }
        }

        avg_order <- new_avg_order
      }
    }
    D[,] <- vapply(pcurves, function(p){ p$dist_ind }, rep(0,nrow(dat)))

    dist_new <- sum(D[W>0], na.rm=TRUE)
    hasConverged <-
    iter <- iter + 1
  }

  names(pcurves) <- paste('curve',seq_along(pcurves),sep='')

  pcurves
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

.clean_curve <- function(pcurve, W_vec, sample_idx){
  stopifnot(class(pcurve) == "principal_curve")

  # force non-negative distances
  pcurve$dist_ind <- abs(pcurve$dist_ind)
  # force pseudotime to start at 0
  pcurve$lambda <- pcurve$lambda - min(pcurve$lambda, na.rm=TRUE)
  pcurve$w <- W_vec
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
    idx <- which(W[,lin] > 0)

    sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)

    pcurve <- princurve::project_to_curve(dat[idx, ,drop = FALSE], s = curve_list[[lin]])
      # note: princurve::project_to_curve changes the input s
    pcurve$s <- pcurve$s[pcurve$ord, ,drop=FALSE]
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
#' @param data a \code{n} by \code{d} matrix. Here, \code{n} could be a subset of the
#' original dataset
#' @param b kernel bandwith
#'
#' @return a smoothed \code{n} by \code{d} matrix
.smoother_func <- function(lambda, dat, idx = NA, b = 1){
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

  avg_curve <- princurve::project_to_curve(dat, avg)
  W_vec <- rowMeans(sapply(pcurve_list, function(p){ p$w }))
  sample_idx <- sort(unique(unlist(lapply(pcurve_list, function(p){ p$idx}))))
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
  lambda <- pcurve$lambda

  dens <- stats::density(0, bw = 1, kernel = "cosine")
  surv <- list(x = dens$x, y = (sum(dens$y) - cumsum(dens$y))/sum(dens$y))
  box_vals <- graphics::boxplot(lambda[common_idx], plot = FALSE)$stats
  surv$x <- .scale_vector(surv$x, lower = box_vals[1], upper = box_vals[5])
  if(box_vals[1] == box_vals[5]){
    pct_l <- rep(0, length(lambda))
  } else {
    pct_l <- stats::approx(surv$x, surv$y, lambda, rule = 2)$y
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

.check_shrinkage <- function(pct_shrink){
  # check for degenerate case (if one curve won't be
  # shrunk, then the other curve shouldn't be,
  # either)
  all_zero <- vapply(pct_shrink[[i]], function(pij){
    return(all(pij == 0))
  }, TRUE)
  if(any(all_zero)){
    if(allow_breaks){
      new_avg_order[[i]] <- NULL
      message('Curves for ', ns[1], ' and ',
              ns[2], ' appear to be going in opposite ',
              'directions. No longer forcing them to ',
              'share an initial point. To manually ',
              'override this, set allow.breaks = ',
              'FALSE.')
    }
    pct_shrink[[i]] <- lapply(pct_shrink[[i]],
                              function(pij){
                                pij[] <- 0
                                return(pij)
                              })
  }
}

.shrink_to_avg <- function(pcurve, avg_curve, pct, dat, stretch = 2){
  n <- nrow(pcurve$s)
  p <- ncol(pcurve$s)
  lam <- pcurve$lambda
  s <- vapply(seq_len(p),function(jj){
    orig_jj <- pcurve$s[,jj]
    avg_jj <- approx(x = avg_curve$lambda, y = avg_curve$s[,jj], xout = lam,
                     rule = 2)$y
    avg_jj * pct + orig_jj * (1-pct)
  }, rep(0,n))
  w <- pcurve$w
  pcurve <- princurve::project_to_curve(dat, as.matrix(s[pcurve$ord, ,drop = FALSE]),
                                        stretch = stretch)
  pcurve$w <- w

  pcurve
}
