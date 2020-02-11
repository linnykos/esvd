.hist_augment <- function(x, breaks = 50, max_val = NA, min_val = 0,
                          lwd = 1, ...){
  if(is.na(max_val)) max_val <- max(x) else x <- x[x <= max_val]
  break_vec <- seq(min_val, max_val, length.out = breaks)

  min_nonzero <- min(x[x != min_val])
  interval <- diff(break_vec)[1]

  if(interval + min_val > min_nonzero){
    interval <- min_nonzero - min_val
  }

  break_vec <- break_vec + interval/2
  break_vec <- c(break_vec[1] - diff(break_vec)[1], break_vec)

  zz <- graphics::hist(x, breaks = break_vec, col = "gray", ...)
  len <- length(which(x == min_val))
  if(len != 0){
    graphics::rect(zz$breaks[1], 0, zz$breaks[2], len, col = grDevices::rgb(0.803, 0.156, 0.211, 0.5))
  }

  invisible()
}
