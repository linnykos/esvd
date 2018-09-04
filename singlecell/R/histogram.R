.hist_augment <- function(x, breaks = 50, min_val = 0,
                         multiplier = 1, param_list = NA, lwd = 1, ...){
  break_vec <- seq(min_val, max(x), length.out = breaks)

  min_nonzero <- min(x[x != min_val])
  interval <- diff(break_vec)[1]

  if(interval + min_val > min_nonzero){
    interval <- min_nonzero - min_val
  }

  break_vec <- break_vec + interval/2
  break_vec <- c(break_vec[1] - diff(break_vec)[1], break_vec)

  zz <- hist(x, breaks = break_vec, col = "gray", ...)
  len <- length(which(x == min_val))
  if(len != 0){
    rect(zz$breaks[1], 0, zz$breaks[2], len, col = rgb(0.803, 0.156, 0.211, 0.5))
  }

  if(!all(is.na(param_list))){
    x_seq <- seq(min_val, max(x), length.out = 1000)

    density_list <- lapply(param_list, function(y){
      lapply(1:2, function(z){ multiplier*likelihood(y[[z]], x_seq) })
    })

    ymax <- max(unlist(density_list), na.rm = T)

    for(i in 1:length(density_list)){
      lines(x_seq, density_list[[i]][[1]], col = rgb(0.803, 0.156, 0.211), lwd = lwd, lty = 1, ... )
      lines(x_seq, density_list[[i]][[2]], col = rgb(0.584, 0.858, 0.564), lwd = lwd, lty = 1, ... )
    }
  }

  invisible()
}
