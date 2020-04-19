#' Plotting diagnostic to determine goodness of fit
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param nat_mat_list list of natural parameter matrices, each of same dimension as \code{dat}
#' @param family character such as \code{"gaussian"}, \code{"exponential"}, \code{"neg_binom"} or \code{"curved_gaussian"}
#' @param missing_idx_list list of missing indices, same length as \code{nat_mat_list}
#' @param width parameter, controlling quantile of prediction region
#' @param scalar additional parameter needed to compute distribution corresponding to \code{family}
#' @param plot boolean
#' @param max_points maximum number of points to be shown in the scatterplot, purely for visualization purposes only
#' @param tol parameter between \code{0} and \code{1} for how strict (\code{1} being the strictest) to measure
#' if the principal angle falls within the prediction region
#' @param ... additional plotting parameters
#'
#' @return either nothing if \code{plot} is \code{TRUE} (and a plot is shown) or the principle angle otherwise
#' @export
plot_prediction_against_observed <- function(dat, nat_mat_list, family, missing_idx_list = list(1:prod(dim(dat))),
                                             width = 0.8, scalar = NA, plot = T,
                                             max_points = 500000, tol = 0.95, ...){
  stopifnot(length(nat_mat_list) == length(missing_idx_list))

  pred_mat_list <- lapply(nat_mat_list, function(nat_mat){
    compute_mean(nat_mat, family = family, scalar = scalar)
  })

  tmp_list <- lapply(1:length(nat_mat_list), function(i){
    cbind(dat[missing_idx_list[[i]]], pred_mat_list[[i]][missing_idx_list[[i]]])
  })

  # compute the principal angle and
  angle_vec <- sapply(tmp_list, .compute_principal_angle)
  angle_val <- mean(angle_vec)

  # determine if the principal angle falls within the prediction region
  tmp_mat <- do.call(rbind, tmp_list)
  colnames(tmp_mat) <- c("observed_val", "predicted_val")
  res <- .within_prediction_region(max(tmp_mat[,"predicted_val"]), family = family, width = width,
                                   scalar = scalar, angle_val = angle_val, tol = tol)

  if(nrow(tmp_mat) > max_points){
    tmp_mat <- tmp_mat[sample(1:nrow(tmp_mat), max_points),]
  }

  if(plot){
    .plot_pca_diagnostic(tmp_mat, seq_vec = res$seq_vec, interval_mat = res$interval_mat,
                         principal_line = res$principal_line, angle_val = angle_val, ...)
  } else {
    list(angle_val = angle_val, bool = res$bool)
  }
}

#' Select tuning parameter for a bunch of matrix completed fits
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param nat_mat_list_list list of lists of natural parameter matrices, each of same dimension as \code{dat}.
#' The length of the list \code{nat_mat_list_list} is equal to \code{\length{scalar_vec}},
#' the number of different scalars to select from, and is
#' comprised of individual lists. The length of each of such lists within each element of
#' \code{nat_mat_list_list} is equal to the length of \code{length(missing_idx_list)}, the
#' number of different trials for each parameter setting
#' @param family character such as \code{"gaussian"}, \code{"exponential"}, \code{"neg_binom"} or \code{"curved_gaussian"}
#' @param missing_idx_list list of missing indices, same length as \code{nat_mat_list_list[[1]]}
#' @param width parameter, controlling quantile of prediction region
#' @param scalar_vec vector of additional parameters needed to compute distribution corresponding to \code{family},
#' of length equal to \code{length(nat_mat_list_list)}
#'
#' @return a list
#' @export
tuning_select_scalar <- function(dat, nat_mat_list_list, family, missing_idx_list = list(1:prod(dim(dat))),
                          width = 0.8, scalar_vec = rep(NA, length(nat_mat_list_list))){
  stopifnot(length(nat_mat_list_list) == length(scalar_vec))
  stopifnot(length(unique(sapply(nat_mat_list_list, length))) == 1)
  stopifnot(length(nat_mat_list_list[[1]]) == length(missing_idx_list))

  res_list <- lapply(1:length(nat_mat_list_list), function(i){
    plot_prediction_against_observed(dat, nat_mat_list_list[[i]], family = family,
                                     missing_idx_list = missing_idx_list,
                                     width = width, scalar = scalar_vec[i], plot = F)
  })

  bool_vec <- sapply(res_list, function(x){x$bool})

  if(any(bool_vec)){
    quality_vec <- sapply(res_list[which(bool_vec)], function(x){x$angle_val})
    scalar_vec2 <- scalar_vec[which(bool_vec)]
  } else {
    quality_vec <- sapply(res_list, function(x){x$angle_val})
    scalar_vec2 <- scalar_vec
  }

  idx <- which.min(abs(quality_vec - 45))
  all_results <- cbind(sapply(res_list, function(x){x$angle_val}),
                       sapply(res_list, function(x){x$bool}),
                       scalar_vec)
  colnames(all_results)
  list(scalar = scalar_vec2[idx], quality = quality_vec[idx], idx = idx, all_results = all_results)
}

#########

.compute_principal_angle <- function(tmp_mat){
  pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
  vec <- pca_res$rotation[,1]; vec <- vec/.l2norm(vec)
  if(sign(vec[1]) < 0)  vec <- -1*vec
  angle_val <- as.numeric(acos(as.numeric(c(0,1) %*% vec)))
  angle_val * 180/pi
}

.within_prediction_region <- function(max_val, family, width, scalar, angle_val, tol = 0.95){
  seq_vec <- seq(0, max_val, length.out = 500)

  interval_mat <- sapply(seq_vec, function(x){
    .compute_prediction_interval_from_mean(x, family = family, width = width, scalar = scalar)
  })
  rownames(interval_mat) <- c("lower", "upper")

  principal_line <- seq_vec * tan(angle_val*pi/180)

  bool_vec <- apply(cbind(interval_mat["lower",] <= principal_line,
                          interval_mat["upper",] >= principal_line), 1, all)
  bool <- sum(bool_vec)/length(bool_vec) >= tol

  list(seq_vec = seq_vec, interval_mat = interval_mat, principal_line = principal_line, bool = bool)
}

.plot_pca_diagnostic <- function(tmp_mat, seq_vec, interval_mat, principal_line, angle_val, ...){
  stopifnot(ncol(interval_mat) == length(principal_line))
  rad <- 2/5*max(tmp_mat[,"predicted_val"])
  seq_max <- 2*max(tmp_mat)

  graphics::plot(NA, asp = T, xlim = range(tmp_mat), ylim = range(tmp_mat),
                 xlab = "Predicted value", ylab = "Observed value", ...)

  graphics::polygon(c(seq_vec, rev(seq_vec)), c(interval_mat["upper",], rev(interval_mat["lower",])), col = grDevices::rgb(1,0,0,0.2),
                    border = NA, density = 30, angle = -45)
  graphics::points(tmp_mat[,"predicted_val"], tmp_mat[,"observed_val"], pch = 16, col = grDevices::rgb(0,0,0,0.2))

  graphics::lines(rep(0, 2), c(-2*seq_max, 2*seq_max), col = "red", lwd = 1)
  graphics::lines(c(-2*seq_max, 2*seq_max), rep(0, 2), col = "red", lwd = 1)
  graphics::lines(c(-2*seq_max, 2*seq_max), c(-2*seq_max, 2*seq_max), col = "red", lwd = 2)

  graphics::lines(seq_vec, interval_mat["lower",], col = "red", lty = 2, lwd = 2)
  graphics::lines(seq_vec, interval_mat["upper",], col = "red", lty = 2, lwd = 2)

  graphics::lines(seq_vec, principal_line, col = "blue", lwd = 2, lty = 2)

  radian_seq <- seq(0, angle_val*pi/180, length.out = 100)
  x_circ <- rad * cos(radian_seq)
  y_circ <- rad * sin(radian_seq)
  graphics::lines(x_circ, y_circ, lty = 2)
  graphics::text(x = rad , y = 2/5*rad, pos = 4, label = paste0(round(angle_val, 1), " degrees"))

  invisible()
}
