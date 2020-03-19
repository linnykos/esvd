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

goodness_heuristic <- function(dat, nat_mat, family, missing_idx = 1:prod(dim(dat)),
                               seq_length = 21, ...){
  pred_mat <- compute_mean(nat_mat, family = family, ...)

  tmp_mat <- cbind(dat[missing_idx], pred_mat[missing_idx])

  width_seq <- seq(0, 1, length = seq_length)

  vec <- sapply(width_seq, function(width){
    interval_mat <- sapply(tmp_mat[,2], function(x){
      .compute_prediction_interval_from_mean(x, family = family, width = width, ...)
    })

    sum(sapply(1:nrow(tmp_mat), function(x){
      interval_mat[1,x] <= tmp_mat[x,1] & interval_mat[2,x] >= tmp_mat[x,1]
    }))/nrow(tmp_mat)
  })
}

plot_prediction_against_observed <- function(dat, nat_mat, family, missing_idx = 1:prod(dim(dat)),
                                             seq_max = NA, width = 0.8, scalar = NA, ...){
  pred_mat <- compute_mean(nat_mat, family = family, scalar = scalar)
  tmp_mat <- cbind(dat[missing_idx], pred_mat[missing_idx])

  if(is.na(seq_max)) seq_max <- max(c(pred_mat, dat))
  seq_vec <- seq(0, seq_max, length.out = 500)

  interval_mat <- sapply(seq_vec, function(x){
    .compute_prediction_interval_from_mean(x, family = family, width = width, scalar = scalar)
  })

  graphics::plot(NA, asp = T, xlim = range(tmp_mat), ylim = range(tmp_mat),
       xlab = "Predicted value", ylab = "Observed value", ...)

  graphics::polygon(c(seq_vec, rev(seq_vec)), c(interval_mat[2,], rev(interval_mat[1,])), col = rgb(1,0,0,0.2),
          border = NA, density = 30, angle = -45)
  graphics::points(tmp_mat[,2], tmp_mat[,1], pch = 16, col = rgb(0,0,0,0.2))

  graphics::lines(rep(0, 2), c(-2*seq_max, 2*seq_max), col = "red", lwd = 1)
  graphics::lines(c(-2*seq_max, 2*seq_max), rep(0, 2), col = "red", lwd = 1)
  graphics::lines(c(-2*seq_max, 2*seq_max), c(-2*seq_max, 2*seq_max), col = "red", lwd = 2)

  graphics::lines(seq_vec, interval_mat[1,], col = "red", lty = 2, lwd = 2)
  graphics::lines(seq_vec, interval_mat[2,], col = "red", lty = 2, lwd = 2)

  pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
  graphics::lines(c(0, 2*seq_max), c(0, 2*seq_max*pca_res$rotation[1,1]/pca_res$rotation[2,1]), col = "blue", lwd = 2, lty = 2)

  rad <- 2/5*max(tmp_mat[,1])
  ang <- as.numeric(acos(abs(c(0,1) %*% pca_res$rotation[,1])))
  radian_seq <- seq(0, ang, length.out = 100)
  x_circ <- rad * cos(radian_seq)
  y_circ <- rad * sin(radian_seq)
  graphics::lines(x_circ, y_circ, lty = 2)
  graphics::text(x = rad , y = 2/5*rad, pos = 4, label = paste0(round(ang * 180/pi, 1), " degrees"))

  invisible()
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

.compute_variance <- function(mean_mat, family, scalar, ...){
  if(family %in% c("neg_binom", "poisson")){
    mean_mat + mean_mat^2/scalar
  } else {
    # this means it's curved gaussian
    (mean_mat/scalar)^2
  }
}
