#' Plotting diagnostic to determine goodness of fit
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param nat_mat_list list of natural parameter matrices, each of same dimension as \code{dat}
#' @param family character such as \code{"gaussian"}, \code{"exponential"}, \code{"neg_binom"} or \code{"curved_gaussian"}
#' @param missing_idx_list list of missing indices, same length as \code{nat_mat_list}
#' @param width parameter, controlling quantile of prediction region. The default is \code{0.8}, meaning
#' that the prediction region is by default from the 10th to 90th quantile.
#' @param scalar additional parameter needed to compute distribution corresponding to \code{family}
#' @param plot boolean
#' @param compute_percentage boolean to whether or not compute the percentage of points that fall within the \code{width}-sized
#' prediction region. It is suggested to keep this as \code{FALSE} since this might be slow to compute.
#' @param max_points maximum number of points to be shown in the scatterplot, purely for visualization purposes only
#' @param tol parameter between \code{0} and \code{1} for how strict (\code{1} being the strictest) to measure
#' if the principal angle falls within the prediction region
#' @param xlim plotting parameter
#' @param ylim plotting parameter
#' @param transparency plotting parameter
#' @param cex_text plotting parameter
#' @param ... additional plotting parameters
#'
#' @return either nothing if \code{plot} is \code{TRUE} (and a plot is shown) or a
#' list of elements include \code{angle_val} (the average principal angle among the \code{length(nat_mat_list)} estimates),
#' \code{angle_sd} (the standard deviation of the principal angle among the \code{length(nat_mat_list)} estimates),
#' \code{bool} (whether or not the average principal angle in \code{angle_val} falls within the
#' \code{width}-sized prediction region) and \code{percentage} (the averge percentage of values
#' that fall within the \code{width}-sized prediction region)
#' @export
plot_prediction_against_observed <- function(dat, nat_mat_list, family, missing_idx_list = list(1:prod(dim(dat))),
                                             width = 0.8, scalar = NA, plot = T,
                                             compute_percentage = F,
                                             max_points = 500000, tol = 0.95, xlim = NA,
                                             ylim = NA, transparency = 0.2, cex_text = 1, ...){
  stopifnot(length(nat_mat_list) == length(missing_idx_list))

  nat_mat_list <- lapply(nat_mat_list, function(nat_mat){
    compute_mean(nat_mat, family = family, scalar = scalar)
  })

  tmp_list <- lapply(1:length(nat_mat_list), function(i){
    cbind(dat[missing_idx_list[[i]]], nat_mat_list[[i]][missing_idx_list[[i]]])
  })

  # compute percentage of points that fall within the region
  if(compute_percentage){
    percentage_vec <- sapply(1:length(nat_mat_list), function(i){
      .compute_overlap_percentage(dat[missing_idx_list[[i]]], nat_mat_list[[i]][missing_idx_list[[i]]],
                                  family = family, width = width, scalar = scalar)
    })
    percentage <- mean(percentage_vec)
  } else{
    percentage <- NA
  }

  # compute the principal angle and
  angle_vec <- sapply(tmp_list, compute_principal_angle)
  angle_val <- mean(angle_vec)
  angle_sd <- stats::sd(angle_vec)

  # determine if the principal angle falls within the prediction region
  tmp_mat <- do.call(rbind, tmp_list)
  colnames(tmp_mat) <- c("observed_val", "predicted_val")
  res <- .within_prediction_region(1.5*max(tmp_mat), family = family, width = width,
                                   scalar = scalar, angle_val = angle_val, tol = tol,
                                   effective_max = max(tmp_mat[,"predicted_val"]))

  if(nrow(tmp_mat) > max_points){
    tmp_mat <- tmp_mat[sample(1:nrow(tmp_mat), max_points),]
  }

  if(plot){
    .plot_pca_diagnostic(tmp_mat, seq_vec = res$seq_vec, interval_mat = res$interval_mat,
                         principal_line = res$principal_line, angle_val = angle_val,
                         xlim = xlim, ylim = ylim, transparency = transparency,
                         cex_text = cex_text, ...)
  } else {
    list(angle_val = angle_val, angle_sd = angle_sd, bool = res$bool, percentage = percentage)
  }
}

#' Select tuning parameter for a bunch of matrix completed fits
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param nat_mat_list_list list of lists of natural parameter matrices, each of same dimension as \code{dat}.
#' The length of the list \code{nat_mat_list_list} is equal to \code{length{scalar_vec}},
#' the number of different scalars to select from, and is
#' comprised of individual lists. The length of each of such lists within each element of
#' \code{nat_mat_list_list} is equal to the length of \code{length(missing_idx_list)}, the
#' number of different trials for each parameter setting
#' @param family character such as \code{"gaussian"}, \code{"exponential"}, \code{"neg_binom"} or \code{"curved_gaussian"}
#' @param missing_idx_list list of missing indices, same length as \code{nat_mat_list_list[[1]]}
#' @param width parameter, controlling quantile of prediction region
#' @param scalar_vec vector of additional parameters needed to compute distribution corresponding to \code{family},
#' of length equal to \code{length(nat_mat_list_list)}
#' @param compute_percentage boolean to whether or not compute the percentage of points that fall within the \code{width}-sized
#' prediction region. It is suggested to keep this as \code{FALSE} since this might be slow to compute.
#'
#' @return a list
#' @export
tuning_select_scalar <- function(dat, nat_mat_list_list, family, missing_idx_list = list(1:prod(dim(dat))),
                          width = 0.8, scalar_vec = rep(NA, length(nat_mat_list_list)),
                          compute_percentage = F){
  stopifnot(length(nat_mat_list_list) == length(scalar_vec))
  stopifnot(length(unique(sapply(nat_mat_list_list, length))) == 1)
  stopifnot(length(nat_mat_list_list[[1]]) == length(missing_idx_list))

  training_idx_list <- lapply(1:length(missing_idx_list), function(j){
    c(1:prod(dim(dat)))[-missing_idx_list[[j]]]
  })

  # angles for testing
  res_test_list <- lapply(1:length(nat_mat_list_list), function(i){
    plot_prediction_against_observed(dat, nat_mat_list_list[[i]], family = family,
                                     missing_idx_list = missing_idx_list,
                                     width = width, scalar = scalar_vec[i],
                                     compute_percentage = compute_percentage, plot = F)
  })

  # angles for training
  res_train_list <- lapply(1:length(nat_mat_list_list), function(i){
    plot_prediction_against_observed(dat, nat_mat_list_list[[i]], family = family,
                                     missing_idx_list = training_idx_list,
                                     width = width, scalar = scalar_vec[i],
                                     compute_percentage = compute_percentage, plot = F)
  })

  # compile all the results
  all_results <- cbind(sapply(res_train_list, function(x){x$angle_val}),
                       sapply(res_train_list, function(x){x$bool}),
                       sapply(res_train_list, function(x){x$percentage}),
                       sapply(res_test_list, function(x){x$angle_val}),
                       sapply(res_test_list, function(x){x$bool}),
                       sapply(res_test_list, function(x){x$percentage}),
                       scalar_vec)
  colnames(all_results) <- c("training_angle", "training_bool",
                             "training_percentage",
                             "testing_angle", "testing_bool",
                             "testing_percentage",
                             "scalar")

  # determine the best parameter
  bool_vec <- all_results[,"testing_bool"]

  if(any(bool_vec == 1)){
    quality_vec <- all_results[which(bool_vec == 1), "testing_angle"]
    scalar_vec2 <- scalar_vec[which(bool_vec == 1)]
    idx_vec <- c(1:length(scalar_vec))[which(bool_vec == 1)]
  } else {
    quality_vec <- all_results[, "testing_angle"]
    scalar_vec2 <- scalar_vec
    idx_vec <- c(1:length(scalar_vec))
  }

  idx <- which.min(abs(quality_vec - 45))

  colnames(all_results)
  list(scalar = scalar_vec2[idx], quality = quality_vec[idx], idx = idx_vec[idx], all_results = all_results)
}

#' Compute principal angle
#'
#' @param tmp_mat a matrix with \code{n} rows (for \code{n} samples) and \code{2} columns,
#' where the first column represents the observed data and the second column represents its
#' corresponding predicted values
#'
#' @return numeric
#' @export
compute_principal_angle <- function(tmp_mat){
  pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
  vec <- pca_res$rotation[,1]; vec <- vec/.l2norm(vec)
  if(sign(vec[1]) < 0)  vec <- -1*vec
  angle_val <- as.numeric(acos(as.numeric(c(0,1) %*% vec)))
  angle_val * 180/pi
}

#########

.compute_overlap_percentage <- function(dat_vec, nat_vec, family, width, scalar){
  interval_mat <- sapply(nat_vec, function(x){
    .compute_prediction_interval_from_mean(x, family = family, width = width, scalar = scalar)
  })

  sum(sapply(1:ncol(interval_mat), function(i){
    dat_vec[i] >= interval_mat[1,i] & dat_vec[i] <= interval_mat[2,i]
  }))/ncol(interval_mat)
}

.within_prediction_region <- function(max_val, family, width, scalar, angle_val, tol = 0.95,
                                      effective_max = max_val){
  seq_vec <- seq(0, max_val, length.out = 500)
  stopifnot(any(seq_vec <= effective_max))

  interval_mat <- sapply(seq_vec, function(x){
    .compute_prediction_interval_from_mean(x, family = family, width = width, scalar = scalar)
  })
  rownames(interval_mat) <- c("lower", "upper")

  principal_line <- seq_vec * tan(angle_val*pi/180)
  idx <- which(seq_vec <= effective_max)

  bool_vec <- apply(cbind(interval_mat["lower",idx] <= principal_line[idx],
                          interval_mat["upper",idx] >= principal_line[idx]), 1, all)
  bool <- sum(bool_vec)/length(bool_vec) >= tol

  list(seq_vec = seq_vec, interval_mat = interval_mat, principal_line = principal_line, bool = bool)
}

.plot_pca_diagnostic <- function(tmp_mat, seq_vec, interval_mat, principal_line, angle_val,
                                 xlim = NA, ylim = NA, transparency = 0.2, cex_text = 1, ...){
  stopifnot(ncol(interval_mat) == length(principal_line))

  rad <- 2/5*max(tmp_mat)
  seq_max <- 2*max(tmp_mat)
  lim_vec <- range(c(0,tmp_mat))
  if(all(is.na(xlim))) xlim <- lim_vec
  if(all(is.na(ylim))) ylim <- lim_vec

  graphics::plot(NA, asp = T, xlim = xlim, ylim = ylim,
                 xlab = "Predicted value", ylab = "Observed value", ...)
  graphics::polygon(c(seq_vec, rev(seq_vec)), c(interval_mat["upper",], rev(interval_mat["lower",])), col = grDevices::rgb(1,0,0,0.2),
                    border = NA, density = 30, angle = -45)
  graphics::points(tmp_mat[,"predicted_val"], tmp_mat[,"observed_val"], pch = 16, col = grDevices::rgb(0,0,0,transparency))

  graphics::lines(rep(0, 2), c(-2*seq_max, 2*seq_max), col = "red", lwd = 1)
  graphics::lines(c(-2*seq_max, 2*seq_max), rep(0, 2), col = "red", lwd = 1)
  graphics::lines(c(-2*seq_max, 2*seq_max), c(-2*seq_max, 2*seq_max), col = "red", lwd = 2)

  graphics::lines(seq_vec, interval_mat["lower",], col = "red", lty = 2, lwd = 2)
  graphics::lines(seq_vec, interval_mat["upper",], col = "red", lty = 2, lwd = 2)

  graphics::lines(seq_vec, principal_line, col = "blue", lwd = 2, lty = 2)

  radian_seq <- seq(0, angle_val*pi/180, length.out = 100)
  x_circ <- rad * cos(radian_seq)
  y_circ <- rad * sin(radian_seq)
  graphics::lines(x_circ, y_circ, lwd = 2, col = "white")
  graphics::lines(x_circ, y_circ, lty = 2)
  graphics::text(x = rad , y = 2/5*rad, pos = 4, label = paste0(round(angle_val, 1), " degrees"),
                 cex = cex_text)

  invisible()
}
