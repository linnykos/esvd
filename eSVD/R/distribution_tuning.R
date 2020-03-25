#' Plotting diagnostic to determine goodness of fit
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param nat_mat_list list of natural parameter matrices, each of same dimension as \code{dat}
#' @param family character such as \code{"gaussian"} or \code{"exponential"}
#' @param missing_idx_list list of missing indices, same length as \code{nat_mat_list}
#' @param seq_max plotting parameter, controlling x- and y-axes
#' @param width plotting parameter, controlling quantile
#' @param scalar additional parameter needed to compute distribution corresponding to \code{family}
#' @param plot boolean
#' @param ... additional plotting parameters
#'
#' @return either nothing if \code{plot} is \code{TRUE} (and a plot is shown) or the principle angle otherwise
#' @export
plot_prediction_against_observed <- function(dat, nat_mat_list, family, missing_idx_list = list(1:prod(dim(dat))),
                                             seq_max = NA, width = 0.8, scalar = NA, plot = T,
                                             max_points = 500000, ...){
  stopifnot(length(nat_mat_list) == length(missing_idx_list))

  pred_mat_list <- lapply(nat_mat_list, function(nat_mat){
    compute_mean(nat_mat, family = family, scalar = scalar)
  })

  tmp_mat <- do.call(rbind, lapply(1:length(nat_mat_list), function(i){
    cbind(dat[missing_idx_list[[i]]], pred_mat_list[[i]][missing_idx_list[[i]]])
  }))

  pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
  rad <- 2/5*max(tmp_mat[,1])
  ang <- as.numeric(acos(abs(c(0,1) %*% pca_res$rotation[,1])))

  if(nrow(tmp_mat) > max_points){
    tmp_mat <- tmp_mat[sample(1:nrow(tmp_mat), max_points),]
  }

  if(plot){
    if(is.na(seq_max)) seq_max <- 2*max(tmp_mat)
    .plot_pca_diagnostic(tmp_mat, family = family, width = width, scalar = scalar, seq_max = seq_max,
                         pca_res = pca_res, rad = rad, ang = ang, ...)
  } else {
    return(ang * 180/pi)
  }
}

.plot_pca_diagnostic <- function(tmp_mat, family, width, scalar, seq_max, pca_res, rad, ang, ...){

  seq_vec <- seq(0, seq_max, length.out = 500)

  interval_mat <- sapply(seq_vec, function(x){
    .compute_prediction_interval_from_mean(x, family = family, width = width, scalar = scalar)
  })

  graphics::plot(NA, asp = T, xlim = range(tmp_mat), ylim = range(tmp_mat),
                 xlab = "Predicted value", ylab = "Observed value", ...)

  graphics::polygon(c(seq_vec, rev(seq_vec)), c(interval_mat[2,], rev(interval_mat[1,])), col = grDevices::rgb(1,0,0,0.2),
                    border = NA, density = 30, angle = -45)
  graphics::points(tmp_mat[,2], tmp_mat[,1], pch = 16, col = grDevices::rgb(0,0,0,0.2))

  graphics::lines(rep(0, 2), c(-2*seq_max, 2*seq_max), col = "red", lwd = 1)
  graphics::lines(c(-2*seq_max, 2*seq_max), rep(0, 2), col = "red", lwd = 1)
  graphics::lines(c(-2*seq_max, 2*seq_max), c(-2*seq_max, 2*seq_max), col = "red", lwd = 2)

  graphics::lines(seq_vec, interval_mat[1,], col = "red", lty = 2, lwd = 2)
  graphics::lines(seq_vec, interval_mat[2,], col = "red", lty = 2, lwd = 2)

  graphics::lines(c(0, 2*seq_max), c(0, 2*seq_max*pca_res$rotation[1,1]/pca_res$rotation[2,1]), col = "blue", lwd = 2, lty = 2)

  radian_seq <- seq(0, ang, length.out = 100)
  x_circ <- rad * cos(radian_seq)
  y_circ <- rad * sin(radian_seq)
  graphics::lines(x_circ, y_circ, lty = 2)
  graphics::text(x = rad , y = 2/5*rad, pos = 4, label = paste0(round(ang * 180/pi, 1), " degrees"))

  invisible()
}
