# code adapted from https://github.com/kstreet13/slingshot

.get_curves <- function(dat, cluster_labels, lineages, shrink = 1,
                        extend = 'y', reweight = TRUE, reassign = F,
                        thresh = 0.001, maxit = 15, stretch = 2,
                        shrink_method = 'cosine',
                        allow.breaks = TRUE, b = 1){
  stopifnot(shrink >= 0 & shrink <= 1)

  ### setup
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  num_lineage <- length(lineages)
  cluster_vec <- 1:k
  n <- nrow(dat)
  p <- ncol(dat)
  centers <- .compute_cluster_center(dat, cluster_mat)

  W <- .initialize_weight_matrix(cluster_mat, lineages)

  ### determine curve hierarchy
  avg_order <- .initialize_curve_hierarchy(lineages, cluster_vec)

  ### initial curves are piecewise linear paths through the tree


  ### track distances between curves and data points to determine convergence
  dist_new <- sum(abs(D[W>0]), na.rm=TRUE)

  it <- 0
  hasConverged <- FALSE
  while (!hasConverged && it < maxit){
    it <- it + 1
    dist_old <- dist_new

    if(reweight | reassign){
      ordD <- order(D)
      W_prob <- W/rowSums(W)
      WrnkD <- cumsum(W_prob[ordD]) / sum(W_prob) #ERROR?
      Z <- D
      Z[ordD] <- WrnkD #ERROR?
    }

    if(reweight){
      Z_prime <- 1-Z^2
      Z_prime[W==0] <- NA
      W0 <- W
      W <- Z_prime / matrixStats::rowMaxs(Z_prime, na.rm = TRUE) #rowMins(D) / D
      W[is.nan(W)] <- 1 # handle 0/0
      W[is.na(W)] <- 0
      W[W > 1] <- 1
      W[W < 0] <- 0
      W[W0==0] <- 0
    }

    if(reassign){
      # add if z < .5
      idx <- Z < .5
      W[idx] <- 1 #(rowMins(D) / D)[idx]

      # drop if z > .9 and w < .1
      ridx <- matrixStats::rowMaxs(Z, na.rm = TRUE) > .9 &
        matrixStats::rowMins(W, na.rm = TRUE) < .1
      W0 <- W[ridx, ]
      Z0 <- Z[ridx, ]
      W0[!is.na(Z0) & Z0 > .9 & W0 < .1] <- 0
      W[ridx, ] <- W0
    }

    ### predict each dimension as a function of lambda (pseudotime)
    for(l in seq_len(num_lineage)){
      pcurve <- pcurves[[l]]
      sample_idx <- sort(unique(unlist(lapply(as.numeric(lineages[[l]]), function(x){which(cluster_mat[,x] == 1)}))))
      s <-  .smoother_func_better(pcurve$lambda, dat[sample_idx,], idx = sample_idx, b = b)
      new_pcurve <- princurve::project_to_curve(dat[sample_idx,], s = s, stretch = stretch)
      new_pcurve$dist_ind <- abs(new_pcurve$dist_ind)
      new_pcurve$lambda <- new_pcurve$lambda -
        min(new_pcurve$lambda, na.rm = TRUE)
      new_pcurve$w <- W[,l]
      pcurves[[l]] <- new_pcurve
    }
    D[,] <- vapply(1:length(pcurves), function(p){
      sample_idx <- sort(unique(unlist(lapply(as.numeric(lineages[[p]]), function(x){which(cluster_mat[,x] == 1)}))))
      vec <- rep(NA, nrow(dat))
      vec[sample_idx] <- pcurves[[p]]$dist_ind
      vec
    }, rep(0,nrow(dat)))

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
    hasConverged <- (abs((dist_old -
                            dist_new)) <= thresh * dist_old)
  }


  ## STEP 5
  if(reweight | reassign){
    ordD <- order(D)
    W_prob <- W/rowSums(W)
    WrnkD <- cumsum(W_prob[ordD]) / sum(W_prob)
    Z <- D
    Z[ordD] <- WrnkD
  }

  if(reweight){
    Z_prime <- 1-Z^2
    Z_prime[W==0] <- NA
    W0 <- W
    W <- Z_prime / matrixStats::rowMaxs(Z_prime,na.rm = TRUE) #rowMins(D) / D
    W[is.nan(W)] <- 1 # handle 0/0
    W[is.na(W)] <- 0
    W[W > 1] <- 1
    W[W < 0] <- 0
    W[W0==0] <- 0
  }

  if(reassign){
    # add if z < .5
    idx <- Z < .5
    W[idx] <- 1 #(rowMins(D) / D)[idx]

    # drop if z > .9 and w < .1
    ridx <- matrixStats::rowMaxs(Z, na.rm = TRUE) > .9 &
      matrixStats::rowMins(W, na.rm = TRUE) < .1
    W0 <- W[ridx, ]
    Z0 <- Z[ridx, ]
    W0[!is.na(Z0) & Z0 > .9 & W0 < .1] <- 0
    W[ridx, ] <- W0
  }

  for(l in seq_len(num_lineage)){
    class(pcurves[[l]]) <- 'principal_curve'
    pcurves[[l]]$w <- W[,l]
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

#' Initial fit of the curves
#'
#' This method is heavily dependent on \code{princurve::project_to_curve()}.
#'
#' @param dat a \code{n} by \code{d} matrix
#' @param lineages list output of \code{.get_lineages()}
#' @param W a weight matrix that is \code{n} by \code{num_lineage} (number of lineages)
#' @param cluster_mat 0-1 matrix output of \code{.construct_cluster_matrix()} that
#' has \code{n} rows and \code{k} column
#' @param centers matrix output of \code{.compute_cluster_center()} that's
#' \code{k} by \code{d}
#'
#' @return a list that contains the \code{num_lineage} curves as \code{principal_curve}
#' and a distance matrix (\code{D}) that contains the squared distance of each point
#' to its repsective lineage curve
.initial_curve_fit <- function(dat, lineages, W, cluster_mat, centers){
  stopifnot(nrow(dat) == nrow(W), ncol(W) == length(lineages),
            ncol(dat) == ncol(centers), nrow(centers) == ncol(cluster_mat))

  n <- nrow(dat); num_lineage <- length(lineages)
  D <- matrix(NA, nrow = n, ncol = num_lineage)
  cluster_vec <- 1:ncol(cluster_mat)
  pcurves <- list()

  for(lin in 1:num_lineage){
    idx <- which(W[,lin] > 0)
    line_initial <- centers[cluster_vec %in% lineages[[lin]], , drop = FALSE]
    line_initial <- line_initial[match(lineages[[lin]], rownames(line_initial)),,
                                 drop = FALSE]

    # initial fit based on just centers
    curve <- princurve::project_to_curve(dat[idx, ,drop = FALSE],
                                         s = line_initial, stretch = 9999) #note: this changes line_initial
    curve$dist_ind <- abs(curve$dist_ind)

    # more refined fit based on samples of those clusters
    sample_idx <- sort(unique(unlist(lapply(as.numeric(lineages[[lin]]), function(x){
      which(cluster_mat[,x] == 1)
    }))))
    pcurve <- princurve::project_to_curve(dat[sample_idx,], s = curve$s[curve$ord, ,drop=FALSE],
                                          stretch=0)

    # force non-negative distances
    pcurve$dist_ind <- abs(pcurve$dist_ind)
    # force pseudotime to start at 0
    pcurve$lambda <- pcurve$lambda - min(pcurve$lambda,
                                         na.rm=TRUE)
    pcurve$w <- W[sample_idx,lin]
    pcurves[[lin]] <- pcurve
    D[sample_idx,lin] <- abs(pcurve$dist_ind)
  }

  list(pcurves = pcurves, D = D)
}


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


.smoother_func_better <- function(lambda, dat2, idx = NA, b = 1){
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
