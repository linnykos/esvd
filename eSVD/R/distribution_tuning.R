#' Tuning negative binomial or curved Gaussian scalar
#'
#' @param dat \code{n} by \code{p} data matrix
#' @param family character (either \code{"neg_binom"} or \code{"curved_gaussian"})
#' @param iter_max positive intger
#' @param search_min positive intger
#' @param search_max positive intger, larger than \code{search_min}
#' @param search_iter positive integer
#' @param search_grid positive integer
#' @param verbose boolean
#' @param ... all remaining parameters needed for either \code{initialization} or \code{matrix_factorization}
#'
#' @return vector of length \code{search_iter}
#' @export
tuning_scalar <- function(dat, family, iter_max = 5, search_min = 1,
                   search_max = 2000, search_iter = 10, search_grid = 10,
                   verbose = F, ...){
  stopifnot(search_max > search_min)

  # fit initial fit
  family_initial <- ifelse(family == "neg_binom", "poisson", "exponential")
  fit <- .tuning_fit(dat, family = family_initial, scalar = NA, ...)
  if(verbose) print("Finished initial fit")

  # determine initial param
  scalar_vec <- rep(NA, iter_max)
  scalar_vec[1] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat, family = family_initial,
                                        search_min = search_min,
                                        search_max = search_max,
                                        search_iter = search_iter,
                                        search_grid = search_grid, ...)
  if(verbose) print("Finished initial search")

  # iterate between fitting and parameter estimation
  for(i in 2:iter_max){
    fit <- .tuning_fit(dat, family = family, scalar = scalar_vec[i-1], ...)
    if(verbose) print(paste0("Finished fit on iteration ", i))

    scalar_vec[i] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat, family = family,
                                          scalar = scalar_vec[i-1],
                                          search_min = search_min,
                                          search_max = search_max,
                                          search_iter = search_iter,
                                          search_grid = search_grid,...)

    if(verbose) print(paste0("Finished search on iteration ", i))
  }

  scalar_vec
}

################

.tuning_fit <- function(dat, family, scalar, ...){
  init <- initialization(dat, family = family, scalar = scalar, ...)
  fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                          family = family, scalar = scalar, ...)
}


###########

.tuning_param_search <- function(dat, u_mat, v_mat, family,
                                 search_min = 1, search_max = 2000, search_iter = 10,
                                 search_grid = 10, ...){
  stopifnot(ncol(u_mat) == ncol(v_mat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat))
  k <- ncol(u_mat); n <- nrow(dat); p <- ncol(dat)
  df_val <- n*p - (n*k + p*k)

  stopifnot(df_val > 0)

  nat_mat <- u_mat %*% t(v_mat)
  mean_mat_tmp <- compute_mean(nat_mat, family, scalar = 1)
  recompute_mean <- family == "neg_binom"

  lo_val <- search_min
  hi_val <- search_max
  min_val <- Inf

  for(iter in 1:search_iter){
    scalar_seq <- seq(lo_val, hi_val, length.out = search_grid)
    obj_seq <- sapply(scalar_seq, function(scalar){
      .compute_tuning_objective(dat, family, nat_mat, mean_mat_tmp, scalar = scalar,
                                recompute_mean = recompute_mean)
    })

    prev_min <- min_val
    min_val <- scalar_seq[which.min(abs(obj_seq - df_val))]

    width <- abs(hi_val - lo_val)
    lo_val <- max(min_val - width/4, search_min)
    hi_val <- min(min_val + width/4, search_max)
  }

  min_val
}

.compute_tuning_objective <- function(dat, family, nat_mat, mean_mat, scalar, recompute_mean){
  mean_mat <- .recompute_mean(nat_mat, mean_mat, scalar = scalar, recompute_mean = recompute_mean)
  var_mat <- .compute_variance(mean_mat, family, scalar = scalar)
  sum((dat-mean_mat)^2/var_mat)
}

.recompute_mean <- function(nat_mat, mean_mat = NA, scalar = NA, recompute_mean = F){
  if(!recompute_mean){
    stopifnot(all(!is.na(mean_mat)))
    return(mean_mat)
  } else {
    # the nat_mat is for the neg_binom family then
    compute_mean(nat_mat, "neg_binom", scalar = scalar)
  }
}

.compute_variance <- function(mean_mat, family, scalar){
  if(family %in% c("neg_binom", "poisson")){
    mean_mat + mean_mat^2/scalar
  } else {
    # this means it's curved gaussian
    (mean_mat/scalar)^2
  }
}
