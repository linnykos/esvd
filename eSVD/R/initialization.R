#' Initialization for matrix factorization
#'
#' This function uses \code{eSVD:::.projected_gradient_descent} primarily, so see
#' that function's documentation for more details.
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param k positive integer less than \code{min(c(nrow(dat), ncol(dat)))}
#' @param family character (\code{"gaussian"}, \code{"exponential"}, \code{"poisson"}, \code{"neg_binom"},
#' or \code{"curved gaussian"})
#' @param max_val maximum magnitude of the inner product
#' @param max_iter numeric
#' @param tol numeric
#' @param verbose boolean
#' @param ... extra arguments, such as nuisance parameters for \code{"neg_binom"}
#' or \code{"curved gaussian"} for \code{family}
#'
#' @return list containing \code{u_mat} and \code{v_mat}
#' @export
initialization <- function(dat, k = 2, family,
                           max_val = NA,
                           max_iter = 10, tol = 1e-3,
                           verbose = F, ...){
  stopifnot(is.matrix(dat), k <= min(dim(dat)), k %% 1 == 0, k > 0)

  direction <- .dictate_direction(family)
  if(!is.na(max_val)){
    stopifnot(max_val > 0)

    # flip max_val so it's the correct sign downstream
    if(!is.na(direction) && direction == "<=") max_val <- -max_val
  }

  # initialize
  dat <- .matrix_completion(dat, k = k)
  attr(dat, "family") <- family

  # projected gradient descent
  nat_mat <- .projected_gradient_descent(dat, k = k, max_val = max_val,
                                          direction = direction,
                                          max_iter = max_iter,
                                          tol = tol, ...)

  res <- .svd_projection(nat_mat, k = k, factors = T)
  u_mat <- res$u_mat; v_mat <- res$v_mat

  if(!is.na(direction)){
    if(direction == "<=") {
      stopifnot(all(nat_mat[which(!is.na(dat))] < 0))
    } else {
      stopifnot(all(nat_mat[which(!is.na(dat))] > 0))
    }
  }

  tmp <- .fix_rank_defficiency_initialization(u_mat, v_mat, direction)
  u_mat <- tmp$u_mat; v_mat <- tmp$v_mat

  tmp <- .reparameterize(u_mat, v_mat)
  u_mat <- tmp$u_mat; v_mat <- tmp$v_mat

  list(u_mat = u_mat, v_mat = v_mat)
}

##################################

#' Fill in missing values
#'
#' Uses \code{softImpute::softImpute} to fill in all the possible missing values.
#' This function enforces all the resulting entries to be non-negative.
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param k positive integer less than \code{min(c(nrow(dat), ncol(dat)))}
#'
#' @return a \code{n} by \code{p} matrix
.matrix_completion <- function(dat, k){
  if(any(is.na(dat))){
    lambda0_val <- softImpute::lambda0(dat)
    res <- softImpute::softImpute(dat, rank.max = k, lambda = min(30, lambda0_val/100))
    diag_mat <- .diag_matrix(res$d[1:k])
    pred_naive <- res$u %*% diag_mat %*% t(res$v)
    dat[which(is.na(dat))] <- pred_naive[which(is.na(dat))]
  }

  pmax(dat, 0)
}

#' Initialize the matrix of natural parameters
#'
#' This function first transforms each entry in \code{dat} according to the inverse function that maps
#' natural parameters to their expectation (according to \code{eSVD:::.mean_transformation}) and then
#' uses \code{eSVD:::.project_rank_feasibility} to get a rank-\code{k} approximation of this matrix
#' that lies within the domain of \code{family}
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes.
#' @param k  positive integer less than \code{min(c(nrow(dat), ncol(dat)))}
#' @param family character (\code{"gaussian"}, \code{"exponential"}, \code{"poisson"}, \code{"neg_binom"},
#' or \code{"curved gaussian"})
#' @param max_val maximum magnitude of the inner product
#' @param tol numeric
#' @param ... extra arguments, such as nuisance parameters for \code{"neg_binom"}
#' or \code{"curved gaussian"} for \code{family}
#'
#' @return \code{n} by \code{p} matrix
.determine_initial_matrix <- function(dat, family, k, max_val = NA, tol = 1e-3, ...){
  dat[which(dat <= tol)] <- tol/2
  nat_mat <- .mean_transformation(dat, family, ...)
  direction <- .dictate_direction(family)

  .project_rank_feasibility(nat_mat, k = k, direction = direction,
                                    max_val = max_val)$matrix
}

#' Projected gradient descent for initial matrix
#'
#' After using \code{eSVD:::.determine_initial_matrix} to determine the initial matrix,
#' uses a projected gradient descent strategy to find a rank-\code{k} initial matrix that minimizes
#' the negative log-likelihood according the distribution specified in \code{family}.
#' See the documentation for \code{eSVD:::.adaptive_gradient_step} for more information.
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes.
#' \code{attr(dat, "family")} needs to encode the \code{family} information, as this is how the function
#' uses generics.
#' @param k  positive integer less than \code{min(c(nrow(dat), ncol(dat)))}
#' @param max_val maximum magnitude of the inner product
#' @param direction character either \code{"<="} or \code{">="} or \code{NA} that dictates the domain of the natural
#' parameter for the family in \code{family}
#' @param max_iter maximum number of iterations
#' @param tol numeric
#' @param ... extra arguments, such as nuisance parameters for \code{"neg_binom"}
#' or \code{"curved gaussian"} for \code{family}
#'
#' @return \code{n} by \code{p} matrix
.projected_gradient_descent <- function(dat, k = 2,
                                        max_val = NA, direction = "<=",
                                        max_iter = 50, tol = 1e-3,
                                        ...){
  nat_mat <- .determine_initial_matrix(dat, k = k, family = attr(dat, "family"), max_val = max_val, ...)
  iter <- 1
  new_obj <- .evaluate_objective_mat(dat, nat_mat, ...)
  old_obj <- Inf

  while(abs(new_obj - old_obj) > tol & iter < max_iter){
    old_obj <- new_obj
    gradient_mat <- .gradient_mat(dat, nat_mat, ...)
    new_mat <- .adaptive_gradient_step(dat, nat_mat, gradient_mat, k = k,
                                       max_val = max_val, direction = direction,
                                       ...)

    new_obj <- .evaluate_objective_mat(dat, new_mat, ...)
    nat_mat <- new_mat
    iter <- iter + 1
  }

  nat_mat
}

#' Do an SVD projection
#'
#' Uses \code{RSpectra::svds} to compute the \code{k} leading singular vectors, but
#' sometimes there are numerical instability issues. In case of crashes, the code
#' then uses the default \code{svd} function.
#'
#' @param mat numeric matrix with \code{n} rows and \code{p} columns
#' @param k positive integer less than \code{min(c(n,p))}.
#' @param factors boolean. If \code{TRUE}, return the factors (i.e., a list containing
#' \code{u_mat} and \code{v_mat}). If \code{FALSE}, return the rank-\code{k} matrix directly.
#' @param u_alone boolean. If \code{TRUE}, then place the singular values on \code{v_mat}
#' @param v_alone boolean. If \code{TRUE}, then place the singular values on \code{u_mat}
#'
#' @return a list or numeric matrix, depending on \code{factors}
.svd_projection <- function(mat, k, factors = F,
                            u_alone = F, v_alone = F){
  stopifnot(min(dim(mat)) >= k, any(!u_alone, !v_alone))

  if(min(dim(mat)) > k+2){
    res <- tryCatch({
      # ask for more singular values than needed to ensure stability
      RSpectra::svds(mat, k = k + 2)
    }, error = function(e){
      svd(mat)
    })
  } else {
    res <- svd(mat)
  }

  diag_mat <- .diag_matrix(res$d[1:k])

  if(factors){
    if(u_alone){
      list(u_mat = res$u[,1:k,drop = F],
           v_mat = res$v[,1:k,drop = F]%*%diag_mat)
    } else if(v_alone) {
      list(u_mat = res$u[,1:k,drop = F]%*%diag_mat,
           v_mat = res$v[,1:k,drop = F])
    } else {
      list(u_mat = res$u[,1:k,drop = F]%*%sqrt(diag_mat),
           v_mat = res$v[,1:k,drop = F]%*%sqrt(diag_mat))
    }
  } else {
    res$u[,1:k,drop = F] %*% diag_mat %*% t(res$v[,1:k,drop = F])
  }
}

#' Adaptive projective gradient descent
#'
#' The projective gradient descent is "adaptive" in the sense
#' that it will find an appropriate step size to ensure that after projection,
#' the objective value descreases. This is a heuristic to simply enable
#' reasonable results, not necessarily theoretically justified or computationally
#' efficient.
#'
#' Be aware of how \code{eSVD:::.project_rank_feasibility} works, as if it cannot
#' successfully find a low-rank representation of the current estimate of the matrix, then
#' it will use \code{eSVD:::.sbm_projection} to do the approximation instead.
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes.
#' \code{attr(dat, "family")} needs to encode the \code{family} information, as this is how the function
#' uses generics.
#' @param nat_mat \code{n} by \code{d} matrix
#' @param gradient_mat \code{n} by \code{d} matrix
#' @param k boolean. If \code{TRUE}, then place the singular values on \code{v_mat}
#' @param max_val  maximum magnitude of the inner product, or \code{NA}
#' @param direction character either \code{"<="} or \code{">="} or \code{NA} that dictates the domain of the natural
#' parameter for the family in \code{family}
#' @param stepsize_init numeric
#' @param stepdown_factor numeric
#' @param max_iter numeric
#' @param ... extra arguments, such as nuisance parameters for \code{"neg_binom"}
#' or \code{"curved gaussian"} for \code{family}
#'
#' @return \code{n} by \code{d} matrix
.adaptive_gradient_step <- function(dat, nat_mat, gradient_mat, k,
                                    max_val = NA, direction = "<=",
                                    stepsize_init = 100, stepdown_factor = 2,
                                    max_iter = 20, ...){
  stepsize <- stepsize_init
  init_obj <- .evaluate_objective_mat(dat, nat_mat, ...)
  iter <- 1

  while(iter > max_iter){
    res <- nat_mat - stepsize*gradient_mat
    new_mat <- .project_rank_feasibility(res, direction = direction,
                               max_val = max_val)$matrix

    if(!any(is.na(new_mat))){
      new_obj <- .evaluate_objective_mat(dat, new_mat, ...)

      if(new_obj < init_obj) return(new_mat)
    }

    stepsize <- stepsize/stepdown_factor
    iter <- iter + 1
  }

  # was not able to project
  nat_mat
}

#' Find low-rank and half-space approximate of a matrix
#'
#' Alternating projection heuristic to find the approximate matrix of \code{mat}
#' that both is rank-\code{k} and has all its entries in the halfspace dictated
#' by \code{direction}. This function alternates
#' between using \code{eSVD:::.svd_projection} and \code{eSVD:::.absolute_threshold}.
#'
#' However, if this function does not converge to a solution after \code{max_iter} iterations,
#' it gives up and uses \code{eSVD:::.sbm_projection} instead, which approximates \code{mat}
#' using a stochastic-block-model variant.
#'
#' @param mat numeric matrix with \code{n} rows and \code{p} columns
#' @param k positive integer less than \code{min(c(n,p))}.
#' @param direction character either \code{"<="} or \code{">="} or \code{NA}
#' @param max_val maximum magnitude of each entry in the approximate
#' @param max_iter numeric
#' @param tol numeric
#'
#' @return \code{n} by \code{d} matrix
.project_rank_feasibility <- function(mat, k, direction, max_val = NA,
                                      max_iter = 20, tol = 1e-6){
  stopifnot(!is.na(max_val) | !is.na(direction))
  if(!is.na(max_val) & !is.na(direction)) stopifnot((direction == "<=" & max_val < 0) | (direction == ">=" & max_val > 0))

  mat_org <- mat
  iter <- 1

  while(iter < max_iter){
    res <- .svd_projection(mat, k = k, factors = T)
    mat <- res$u_mat %*% t(res$v_mat)

    if(is.na(direction)){
      if(all(abs(mat) <= max_val+tol)) return(list(matrix = mat, iter = iter))
    } else if (direction == "<=") {
      if(all(mat < 0) && (is.na(max_val) || all(mat > max_val-tol))) return(list(matrix = mat, iter = iter))
    } else {
      if(all(mat > 0) && (is.na(max_val) || all(mat < max_val+tol))) return(list(matrix = mat, iter = iter))
    }

    mat <- .absolute_threshold(mat, direction, max_val)

    iter <- iter + 1
  }

  # if the alternating projection strategy above failed, use a SBM-projection
  mat <- .absolute_threshold(mat_org, direction, max_val)
  res <- .sbm_projection(mat, k)

  list(matrix = res, iter = NA)
}

#' Hard threshold all the entries in a matrix
#'
#' While simple conceptually, this function is a bit janky-looking in its
#' implementation due to numeric instabilities previously encountered.
#'
#' @param mat numeric matrix with \code{n} rows and \code{p} columns
#' @param direction character either \code{"<="} or \code{">="} or \code{NA}
#' @param max_val maximum magnitude of each entry in the approximate
#' @param tol2 numeric
#'
#' @return \code{n} by \code{d} matrix
.absolute_threshold <- function(mat, direction, max_val = NA, tol2 = 1e-3){
  if(is.na(direction)){
    max_val <- abs(max_val)

    idx <- which(abs(mat) >= max_val)
    val <- mat[idx]
    mat[idx] <- sign(val)*(max_val-tol2)

  } else if (direction == "<=") {
    tol <- -tol2
    if(any(mat < 0)) tol <- min(max(mat[mat < 0]), tol)
    stopifnot(tol < 0)
    mat[mat > -tol] <- tol
    if(!is.na(max_val)) mat[mat < max_val] <- max_val+tol2

  } else {
    tol <- tol2
    if(any(mat > 0)) tol <- max(min(mat[mat > 0]), tol)
    stopifnot(tol > 0)
    mat[mat < tol] <- tol
    if(!is.na(max_val)) mat[mat > max_val] <- max_val-tol2
  }

  mat
}

#' Fix rank defficiency among two matrices
#'
#' Given two matrices, \code{u_mat} and \code{v_mat} (both with the same number of columns),
#' adjust these two matrices so the \code{u_mat \%*\% t(v_mat)} is actually of the desired
#' rank. This function is needed since sometimes upstream, the matrices \code{u_mat} and \code{v_mat}
#' do not actually have the rank equal to the number of columns (i.e., empirically we have
#' observed that \code{u_mat} might have a column that is all constant).
#'
#' @param u_mat a numeric matrix
#' @param v_mat a numeric matrix with the same number of columns as \code{u_mat}
#' @param direction character either \code{"<="} or \code{">="} or \code{NA}
#'
#' @return a list of \code{u_mat} and \code{v_mat}
.fix_rank_defficiency_initialization <- function(u_mat, v_mat, direction){
  k <- ncol(u_mat)
  nat_mat <- u_mat %*% t(v_mat)
  k2 <- as.numeric(Matrix::rankMatrix(nat_mat))

  if(k != k2){
    stopifnot(k2 < k)
    sign_val <- ifelse(direction == ">=", 1, -1)

    sd_val <- mean(c(apply(u_mat[,1:k2, drop = F], 2, stats::sd),apply(v_mat[,1:k2, drop = F], 2, stats::sd)))
    for(i in (k2+1):k){
      u_mat[,i] <- abs(stats::rnorm(nrow(u_mat), sd = sd_val))
      v_mat[,i] <- sign_val*abs(stats::rnorm(nrow(v_mat), sd = sd_val))
    }
  }

  list(u_mat = u_mat, v_mat = v_mat)
}
