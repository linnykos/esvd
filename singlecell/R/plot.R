.plot_singlecell <- function(dat, ...){
  col_vec <- .colorRamp_custom(c(0.584, 0.858, 0.564), c(0.803, 0.156, 0.211), 19)
  col_vec <- c("white", col_vec)

  tmp <- as.numeric(dat)
  tmp <- tmp[tmp!=0]

  break_vec <- stats::quantile(tmp, probs = seq(0, 1, length.out = 20))
  break_vec <- c(-5, break_vec)

  graphics::image(.rotate(dat), breaks = break_vec, col = col_vec, asp = nrow(dat)/ncol(dat),
        axes = F, ...)

  invisible()
}


.colorRamp_custom <- function(vec1, vec2, length){
  mat <- matrix(0, nrow = length, ncol = 3)
  for(i in 1:3){
    mat[,i] <- seq(vec1[i], vec2[i], length.out = length)
  }

  luminosity_vec <- apply(mat, 1, function(x){
    0.2126*x[1] + 0.7152*x[2] + 0.0722*x[3]
  })

  target_luminosity <- mean(c(luminosity_vec[1], luminosity_vec[length]))

  mat <- t(sapply(1:nrow(mat), function(x){
    factor <- min(c(target_luminosity/luminosity_vec[x], 1/mat[x,]))
    mat[x,] * factor
  }))

  apply(mat, 1, function(x){
    grDevices::rgb(x[1], x[2], x[3])
  })
}

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

