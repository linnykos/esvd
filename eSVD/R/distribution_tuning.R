#' Tuning negative binomial or curved Gaussian scalar
#'
#' @param dat \code{n} by \code{p} data matrix
#' @param family character (either \code{"neg_binom"} or \code{"curved_gaussian"})
#' @param iter_max positive intger
#' @param search_min positive intger
#' @param search_max positive intger, larger than \code{search_min}
#' @param verbose boolean
#' @param ... all remaining parameters needed for either \code{initialization} or \code{matrix_factorization}
#'
#' @return vector of length \code{search_iter}
#' @export
tuning_scalar <- function(dat, family, iter_max = 5,
                          search_min = 1,
                          search_max = 2000,
                          verbose = F, ...){
  stopifnot(search_max > search_min)
  missing_idx <- eSVD::construct_missing_values(n = nrow(dat), p = ncol(dat))

  # fit initial fit
  family_initial <- ifelse(family == "neg_binom", "poisson", "exponential")
  dat_NA <- dat; dat_NA[missing_idx] <- NA
  missing_val <- dat[missing_idx]
  fit <- .tuning_fit(dat_NA, family = family_initial, scalar = NA, ...)
  if(verbose) print("Finished initial fit")

  # determine initial param
  scalar_vec <- rep(NA, iter_max)
  scalar_vec[1] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat,
                                        family = family_initial,
                                        idx = missing_idx,
                                        search_min = search_min,
                                        search_max = search_max, ...)
  if(verbose) print("Finished initial search")

  # iterate between fitting and parameter estimation
  for(i in 2:iter_max){
    fit <- .tuning_fit(dat_NA, family = family, scalar = scalar_vec[i-1], ...)
    if(verbose) print(paste0("Finished fit on iteration ", i))

    scalar_vec[i] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat, family = family,
                                          idx = missing_idx,
                                          scalar = scalar_vec[i-1],
                                          search_min = search_min,
                                          search_max = search_max, ...)

    if(verbose) print(paste0("Finished search on iteration ", i))
  }

  scalar_vec
}

goodness_heuristic <- function(dat, nat_mat, family, missing_idx = 1:prod(dim(dat)), ...){
  pred_mat <- compute_mean(nat_mat, family = family)

  tmp_mat <- cbind(dat[missing_idx], pred_mat[missing_idx])

  interval_mat <- sapply(tmp_mat[,2], function(x){
    .compute_prediction_interval_from_mean(x, family = family, ...)
  })
}

################

.tuning_fit <- function(dat, family, scalar, ...){
  init <- initialization(dat, family = family, scalar = scalar, ...)
  fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                          family = family, scalar = scalar, ...)
}


###########

.tuning_param_search <- function(dat, u_mat, v_mat, family,
                                 idx = 1:prod(dim(dat)),
                                 search_min = 1, search_max = 2000, ...){
  stopifnot(ncol(u_mat) == ncol(v_mat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat))
  k <- ncol(u_mat); n <- nrow(dat); p <- ncol(dat)
  if(length(idx) == prod(dim(dat))){
    df_val <- n*p - (n*k + p*k)
  } else {
    df_val <- length(idx)
  }

  stopifnot(df_val > 0)

  nat_mat <- u_mat %*% t(v_mat)

  fn <- function(x){
    mean_mat <- compute_mean(nat_mat, family, scalar = x)
    var_mat <- .compute_variance(mean_mat, family, scalar = x)
    abs(sum((dat[idx]-mean_mat[idx])^2/var_mat[idx]) - df_val)
  }

  res <- stats::optimize(fn, interval = c(search_min, search_max))
  res$minimum
}

.compute_variance <- function(mean_mat, family, scalar){
  if(family %in% c("neg_binom", "poisson")){
    mean_mat + mean_mat^2/scalar
  } else {
    # this means it's curved gaussian
    (mean_mat/scalar)^2
  }
}
