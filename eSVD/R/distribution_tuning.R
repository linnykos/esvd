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
                                        search_max = search_max, ...)
  if(verbose) print("Finished initial search")

  # iterate between fitting and parameter estimation
  for(i in 2:iter_max){
    fit <- .tuning_fit(dat, family = family, scalar = scalar_vec[i-1], ...)
    if(verbose) print(paste0("Finished fit on iteration ", i))

    scalar_vec[i] <- .tuning_param_search(dat, fit$u_mat, fit$v_mat, family = family,
                                          scalar = scalar_vec[i-1],
                                          search_min = search_min,
                                          search_max = search_max, ...)

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
                                 search_min = 1, search_max = 2000, ...){
  stopifnot(ncol(u_mat) == ncol(v_mat), nrow(dat) == nrow(u_mat), ncol(dat) == nrow(v_mat))
  k <- ncol(u_mat); n <- nrow(dat); p <- ncol(dat)
  df_val <- n*p - (n*k + p*k)

  stopifnot(df_val > 0)

  nat_mat <- u_mat %*% t(v_mat)
  mean_mat_base <- compute_mean(nat_mat, family, scalar = 1)

  fn <- function(x){
    mean_mat <- compute_mean(nat_mat, family, scalar = x)
    var_mat <- .compute_variance(mean_mat, family, scalar = x)
    abs(sum((dat-mean_mat)^2/var_mat) - df_val)
  }

  gr <- function(x){
    mean_mat <- compute_mean(nat_mat, family, scalar = x)
    var_mat <- .compute_variance(mean_mat, family, scalar = x)
    .compute_gradient(dat, mean_mat, var_mat, df_val, family, scalar = x)
  }

  res <- stats::optim(par = (search_min+search_max)/2, fn = fn, gr = gr,
                      method = "L-BFGS-B",
                      lower = search_min, upper = search_max)
  res$par
}

.compute_variance <- function(mean_mat, family, scalar){
  if(family %in% c("neg_binom", "poisson")){
    mean_mat + mean_mat^2/scalar
  } else {
    # this means it's curved gaussian
    (mean_mat/scalar)^2
  }
}

.compute_gradient <- function(dat, mean_mat, var_mat, df_val, family, scalar){
  sign_val <- sign(sum((dat-mean_mat)^2/var_mat) - df_val)

  if(family %in% c("neg_binom", "poisson")){
    sign_val*sum((dat - mean_mat)^2/var_mat^2*(mean_mat/scalar)^2)
  } else {
    sign_val*(2*scalar*sum((dat - mean_mat)^2/mean_mat^2))
  }
}
