.plot_singlecell <- function(dat, luminosity = T, ...){
  col_vec <- .colorRamp_custom(c(0.584, 0.858, 0.564), c(0.803, 0.156, 0.211), 19,
                               luminosity = luminosity)
  col_vec <- c("white", col_vec)

  tmp <- as.numeric(dat)
  tmp <- tmp[tmp!=0]

  break_vec <- stats::quantile(tmp, probs = seq(0, 1, length.out = 20))
  break_vec <- c(-5, break_vec)

  graphics::image(.rotate(dat), breaks = break_vec, col = col_vec, asp = nrow(dat)/ncol(dat),
        axes = F, ...)

  invisible()
}


.colorRamp_custom <- function(vec1, vec2, length, luminosity = T){
  mat <- matrix(0, nrow = length, ncol = 3)
  for(i in 1:3){
    mat[,i] <- seq(vec1[i], vec2[i], length.out = length)
  }

  if(luminosity){
    luminosity_vec <- apply(mat, 1, function(x){
      0.2126*x[1] + 0.7152*x[2] + 0.0722*x[3]
    })

    target_luminosity <- mean(c(luminosity_vec[1], luminosity_vec[length]))

    mat <- t(sapply(1:nrow(mat), function(x){
      factor <- min(c(target_luminosity/luminosity_vec[x], 1/mat[x,]))
      mat[x,] * factor
    }))
  }

  apply(mat, 1, function(x){
    grDevices::rgb(x[1], x[2], x[3])
  })
}

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

#' Compute the points of the ellipse on a two-dimensional plot
#'
#' Function is for plotting purposes only
#'
#' @param mean_vec 2-dimensional vector
#' @param cov_mat 2 by 2 PSD matrix
#' @param scale_factor scaling factor to determine the size of ellipse
#'
#' @return matrix
.compute_ellipse_points <- function(mean_vec, cov_mat, scale_factor = 1){
  stopifnot(length(mean_vec) == 2, all(dim(cov_mat) == c(2,2)))
  stopifnot(sum(abs(cov_mat - t(cov_mat))) <= 1e-6)

  eig <- eigen(cov_mat)
  alpha <- atan(eig$vectors[2,1]/eig$vectors[1,1])
  if(alpha < 0) alpha <- alpha + 2*pi

  a <- sqrt(eig$values[1])*scale_factor
  b <- sqrt(eig$values[2])*scale_factor

  theta_grid <- seq(0, 2*pi, length.out = 100)
  ellipse_x <- a*cos(theta_grid)
  ellipse_y <- b*sin(theta_grid)

  R <- matrix(c(cos(alpha), -sin(alpha), sin(alpha), cos(alpha)), 2, 2)
  val <- cbind(ellipse_x, ellipse_y) %*% R

  t(apply(val, 1, function(x){x + mean_vec}))
}




