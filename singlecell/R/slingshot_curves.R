# code adapted from https://github.com/kstreet13/slingshot

.get_curves <- function(dat, cluster_mat, lineages, shrink = TRUE,
                        extend = 'y', reweight = TRUE, reassign = TRUE,
                        thresh = 0.001, maxit = 15, stretch = 2,
                        smoother = 'smooth.spline', shrink_method = 'cosine',
                        allow.breaks = TRUE){
  shrink <- as.numeric(shrink)
  # CHECKS
  if(shrink < 0 | shrink > 1){
    stop("'shrink' parameter must be logical or numeric between",
         "0 and 1")
  }

  smoother_func <- .smoother_slingshot(smoother)

  ### setup
  num_lineage <- length(grep("Lineage", names(lineages))) # number of lineages
  clusters <- colnames(cluster_mat)
  n <- nrow(dat)
  p <- ncol(dat)
  nclus <- length(clusters)
  centers <- .compute_clustercenter(dat, cluster_mat)

  if(p == 1){
    centers <- t(centers)
    rownames(centers) <- clusters
  }

  rownames(centers) <- clusters
  W <- vapply(seq_len(num_lineage), function(l){
    rowSums(cluster_mat[, lineages[[num_lineage]], drop = FALSE])
  }, rep(0,nrow(dat))) # weighting matrix
  rownames(W) <- rownames(dat)
  colnames(W) <- names(lineages)[seq_len(num_lineage)]
  W.orig <- W
  D <- W; D[,] <- NA

  ### determine curve hierarchy
  C <- as.matrix(vapply(lineages[seq_len(num_lineage)], function(lin) {
    vapply(clusters, function(cluster_id) {
      as.numeric(cluster_id %in% lin)
    }, 0)
  }, rep(0,nclus)))
  rownames(C) <- clusters
  segmnts <- unique(C[rowSums(C)>1,,drop = FALSE])
  segmnts <- segmnts[order(rowSums(segmnts),decreasing = FALSE), ,
                     drop = FALSE]
  avg_order <- list() #avg_order will have length equal to number of connected components, i think
  for(i in seq_len(nrow(segmnts))){
    idx <- segmnts[i,] == 1
    avg_order[[i]] <- colnames(segmnts)[idx]
    new_col <- rowMeans(segmnts[,idx, drop = FALSE])
    segmnts <- cbind(segmnts[, !idx, drop = FALSE], new_col)
    colnames(segmnts)[ncol(segmnts)] <- paste('average',i,sep='')
  }

  ### initial curves are piecewise linear paths through the tree
  pcurves <- list()
  for(l in seq_len(num_lineage)){
    idx <- which(W[,l] > 0)
    line_initial <- centers[clusters %in% lineages[[l]], ,
                            drop = FALSE]
    line_initial <- line_initial[match(lineages[[l]],
                                       rownames(line_initial)),  ,
                                 drop = FALSE]
    K <- nrow(line_initial)
    # special case: single-cluster lineage
    if(K == 1){
      pca <- stats::prcomp(X[idx, ,drop = FALSE])
      ctr <- line_initial
      line_initial <- rbind(ctr - 10*pca$sdev[1] *
                              pca$rotation[,1], ctr,
                            ctr + 10*pca$sdev[1] *
                              pca$rotation[,1])
      curve <- princurve::project_to_curve(X[idx, ,drop = FALSE],
                                s = line.initial, stretch = 9999)
      # do this twice because all points should have projections
      # on all lineages, but only those points on the lineage
      # should extend it
      pcurve <- princurve::project_to_curve(X, s = curve$s[curve$ord, ,
                                                drop = FALSE], stretch=0)
      pcurve$dist_ind <- abs(pcurve$dist_ind)
      # ^ force non-negative distances
      pcurve$lambda <- pcurve$lambda - min(pcurve$lambda,
                                           na.rm=TRUE)
      # ^ force pseudotime to start at 0
      pcurve$w <- W[,l]
      pcurves[[l]] <- pcurve
      D[,l] <- abs(pcurve$dist_ind)
      next
    }

    if(extend == 'y'){
      curve <- princurve::project_to_curve(dat[idx, ,drop = FALSE],
                                s = line_initial, stretch = 9999) #note: this changes line_initial
      curve$dist_ind <- abs(curve$dist_ind)
    }

    pcurve <- princurve::project_to_curve(dat, s = curve$s[curve$ord, ,drop=FALSE],
                               stretch=0)
    # force non-negative distances
    pcurve$dist_ind <- abs(pcurve$dist_ind)
    # force pseudotime to start at 0
    pcurve$lambda <- pcurve$lambda - min(pcurve$lambda,
                                         na.rm=TRUE)
    pcurve$w <- W[,l]
    pcurves[[l]] <- pcurve
    D[,l] <- abs(pcurve$dist_ind)
  }

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
      s <- pcurve$s
      ordL <- order(pcurve$lambda)
      for(jj in seq_len(p)){
        s[, jj] <- smoother_func(pcurve$lambda, dat[,jj], w = pcurve$w)[ordL]
      }
      new_pcurve <- princurve::project_to_curve(dat, s = s, stretch = stretch)
      new_pcurve$dist_ind <- abs(new_pcurve$dist_ind)
      new_pcurve$lambda <- new_pcurve$lambda -
        min(new_pcurve$lambda, na.rm = TRUE)
      new_pcurve$w <- W[,l]
      pcurves[[l]] <- new_pcurve
    }
    D[,] <- vapply(pcurves, function(p){ p$dist_ind }, rep(0,nrow(dat)))

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

          avg <- .avg_curves(to_avg, dat, stretch = stretch)
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

# DEFINE SMOOTHER FUNCTION
.smoother_slingshot <- function(smoother){
  switch(smoother, loess = function(lambda, xj,
                                    w = NULL, ...){
    loess(xj ~ lambda, weights = w, ...)$fitted
  }, smooth.spline = function(lambda, xj, w = NULL, ..., df = 5,
                              tol = 1e-4){
    # fit <- smooth.spline(lambda, xj, w = w, ..., df = df,
    #                      tol = tol, keep.data = FALSE)
    fit <- tryCatch({
      smooth.spline(lambda, xj, w = w, ..., df = df,
                    tol = tol, keep.data = FALSE)
    }, error = function(e){
      smooth.spline(lambda, xj, w = w, ..., df = df,
                    tol = tol, keep.data = FALSE, spar = 1)
    })
    predict(fit, x = lambda)$y
  })
}

.avg_curves <- function(pcurves, dat, stretch = 2){
  n <- nrow(pcurves[[1]]$s)
  p <- ncol(pcurves[[1]]$s)
  lambdas_all <- lapply(pcurves, function(pcv){pcv$lambda})
  lambdas_all <- unique(unlist(lambdas_all))
  max_shared_lambda <- min(vapply(pcurves, function(pcv){max(pcv$lambda)},0))
  lambdas_all <- sort(lambdas_all[lambdas_all <= max_shared_lambda])

  pcurves_dense <- lapply(pcurves, function(pcv){
    vapply(seq_len(p),function(jj){
      interpolated <- stats::approx(pcv$lambda, pcv$s[,jj], xout = lambdas_all)$y

      interpolated
    }, rep(0,length(lambdas_all)))
  })

  avg <- vapply(seq_len(p),function(jj){
    dim_all <- vapply(seq_along(pcurves_dense),function(i){
      pcurves_dense[[i]][,jj]
    }, rep(0,length(lambdas_all)))

    rowMeans(dim_all)
  }, rep(0,length(lambdas_all)))

  avg_curve <- princurve::project_to_curve(dat, avg, stretch=stretch)
  avg_curve$w <- rowMeans(vapply(pcurves, function(p){ p$w }, rep(0,n)))

  avg_curve
}

.percent_shrinkage <- function(crv, share_idx, method = 'cosine'){
  pst <- crv$lambda
  if(method %in% base::eval(formals(stats::density.default)$kernel)){
    dens <- stats::density(0, bw=1, kernel = method)
    surv <- list(x = dens$x, y = (sum(dens$y) - cumsum(dens$y))/sum(dens$y))
    box_vals <- graphics::boxplot(pst[share_idx], plot = FALSE)$stats
    surv$x <- .scaleAB(surv$x, a = box_vals[1], b = box_vals[5])
    if(box_vals[1]==box_vals[5]){
      pct_l <- rep(0, length(pst))
    }else{
      pct_l <- stats::approx(surv$x, surv$y, pst, rule = 2)$y
    }
  }

  pct_l
}

.scaleAB <- function(x,a=0,b=1){
  ((x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))*(b-a)+a
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

.smoother_func_better <- function(lambda, dat2){
  ord <- order(lambda, decreasing = F)
  lambda <- lambda[ord]
  dat2 <- dat2[ord,]

  kernel_func <- function(x, y){exp(-(x-y)^2)}

  t(sapply(1:length(lambda), function(i){
    weights <- kernel_func(lambda[i], lambda)
    as.numeric(colSums(diag(weights) %*% dat2))/sum(weights)
  }))
}
