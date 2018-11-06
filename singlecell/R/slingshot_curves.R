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
#' @param shrink_method
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
                        shrink_method = 'cosine',
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
      .smoother_func(pcurve_list[[lin]]$lambda, dat[sample_idx,],
                     idx = sample_idx, b = b)
    })

    res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)
    pcurve_list <- res$pcurve_list; D <- res$D

    # shrink together lineages near shared clusters
    if(shrink > 0){
      if(max(rowSums(C)) > 1){

        segmnts <- unique(C[rowSums(C)>1,,drop=FALSE])
        segmnts <- segmnts[order(rowSums(segmnts), decreasing = FALSE),, drop = FALSE]
        seg_mix <- segmnts
        avg_lines <- list()
        pct_shrink <- list()

        ### determine average curves and amount of shrinkage
        for(i in seq_along(avg_order)){
          ns <- avg_order[[i]]
          to_avg <- lapply(ns, function(n_element){
            if(grepl('Lineage', n_element)){
              l_ind <- as.numeric(gsub('Lineage','',n_element))
              pcurves[[l_ind]]
            } else if(grepl('average', n_element)){
              a_ind <- as.numeric(gsub('average','',n_element))
              avg_lines[[a_ind]]
            }
          })

          avg <- .avg_curves_tmp(to_avg, dat, stretch = stretch)
          avg_lines[[i]] <- avg
          common_ind <- rowMeans(vapply(to_avg, function(crv){ crv$w > 0 },
                                        rep(TRUE,nrow(dat)))) == 1
          pct_shrink[[i]] <- lapply(to_avg, function(crv){
            .percent_shrinkage(crv, common_ind, method = shrink_method)
          })

          # check for degenerate case (if one curve won't be
          # shrunk, then the other curve shouldn't be,
          # either)
          new_avg_order <- avg_order
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

.clean_curve <- function(pcurve, W){
  stopifnot(class(pcurve) == "principal_curve")

  # force non-negative distances
  pcurve$dist_ind <- abs(pcurve$dist_ind)
  # force pseudotime to start at 0
  pcurve$lambda <- pcurve$lambda - min(pcurve$lambda, na.rm=TRUE)
  pcurve$w <- W

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

  pcurve_list <- list()
  for(lin in 1:num_lineage){
    idx <- which(W[,lin] > 0)

    sample_idx <- .determine_idx_lineage(lineages[[lin]], cluster_mat)

    pcurve <- princurve::project_to_curve(dat[idx, ,drop = FALSE], s = curve_list[[lin]])
      # note: princurve::project_to_curve changes the input s
    pcurve$s <- pcurve$s[pcurve$ord, ,drop=FALSE]
    pcurve <- .clean_curve(pcurve, W[sample_idx, lin])
    pcurve_list[[lin]] <- pcurve

    D[sample_idx,lin] <- abs(pcurve$dist_ind)
  }

  list(pcurve_list = pcurve_list, D = D)
}



##################

.avg_curves_tmp <- function(pcurves, dat, stretch = 2){
  n <- nrow(pcurves[[1]]$s)
  p <- ncol(pcurves[[1]]$s)
  lambdas_all <- lapply(pcurves, function(pcv){pcv$lambda})
  lambdas_all <- unique(unlist(lambdas_all))
  max_shared_lambda <- min(vapply(pcurves, function(pcv){max(pcv$lambda)},0))
  lambdas_all <- sort(lambdas_all[lambdas_all <= max_shared_lambda])

  pcurves_dense <- lapply(pcurves, function(pcv){
    vapply(seq_len(p),function(jj){
      stats::approx(pcv$lambda, pcv$s[,jj], xout = lambdas_all)$y
    }, rep(0,length(lambdas_all)))
  })

  avg <- vapply(seq_len(p),function(jj){
    dim_all <- vapply(seq_along(pcurves_dense),function(i){
      pcurves_dense[[i]][,jj]
    }, rep(0,length(lambdas_all)))

    rowMeans(dim_all)
  }, rep(0,length(lambdas_all)))

  avg_curve <- princurve::project_to_curve(dat, avg, stretch=stretch)
  avg_curve$w <- rowMeans(vapply(pcurves, function(p){ p$w }, rep(0,nrow(dat))))

  avg_curve
}


.smoother_func <- function(lambda, dat2, idx = NA, b = 1){
  if(all(!is.na(idx)) & nrow(dat2) != length(lambda)) lambda <- lambda[idx]
  ord <- order(lambda, decreasing = F)
  lambda <- lambda[ord]
  dat2 <- dat2[ord,]

  kernel_func <- function(x, y){exp(-(x-y)^2/b)}

  t(sapply(1:length(lambda), function(i){
    weights <- kernel_func(lambda[i], lambda)
    as.numeric(colSums(diag(weights) %*% dat2))/sum(weights)
  }))
}
