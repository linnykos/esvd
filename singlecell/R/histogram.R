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

  zz <- graphics::hist(x, breaks = break_vec, col = "gray", ...)
  len <- length(which(x == min_val))
  if(len != 0){
    graphics::rect(zz$breaks[1], 0, zz$breaks[2], len, col = grDevices::rgb(0.803, 0.156, 0.211, 0.5))
  }

  if(!all(is.na(param_list))){
    x_seq <- seq(min_val, max(x), length.out = 5000)
    x_seq <- x_seq[-1]

    density_list <- lapply(param_list, function(y){
      tmp <- vector("list", 2)
      tmp[[1]] <- multiplier * y[["proportion"]] * likelihood(y[[1]], x_seq)
      tmp[[2]] <- multiplier * (1-y[["proportion"]]) * likelihood(y[[2]], x_seq)
      tmp
    })

    ymax <- max(unlist(density_list), na.rm = T)

    for(i in 1:length(density_list)){
      graphics::lines(x_seq, density_list[[i]][[1]], col = grDevices::rgb(0.803, 0.156, 0.211), lwd = lwd, lty = 1, ... )
      graphics::lines(x_seq, density_list[[i]][[2]], col = grDevices::rgb(0.584, 0.858, 0.564), lwd = lwd, lty = 1, ... )
    }

    i <- 1
    idx <- min(which(density_list[[i]][[1]] < density_list[[i]][[2]]))
    if(!is.na(idx)){
      graphics::lines(rep(x_seq[idx], 2), c(1e6,-1e6), lwd = lwd, lty = 2, ...)
    }
  }

  invisible()
}


.compute_dropout <- function(param, x){
  like1 <- param[["proportion"]] * likelihood(param[["class1"]], x)
  like2 <- (1-param[["proportion"]]) * likelihood(param[["class2"]], x)

  estimated <- rep(2, length(x))
  estimated[like1 > like2] <- 1

  estimated
}

.compute_dropout_table <- function(param, data_obj){
  estimated <- .compute_dropout(param, data_obj$x)

  table(data_obj$assignment, estimated)
}

.visualize_table <- function(param, data_obj, ...){
  tab <- .compute_dropout_table(param, data_obj)

  col1 <- c(0.584, 0.858, 0.564)
  col2 <- c(0.803, 0.156, 0.211)

  prob1 <- tab[1,1]/sum(tab[1,])
  prob2 <- tab[2,2]/sum(tab[2,])

  col1 <- sapply(col1, function(x){1 - prob1 * (1 - x)})
  col2 <- sapply(col2, function(x){1 - prob2 * (1 - x)})

  graphics::par(mar = c(1,1,3,0))
  graphics::plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T, axes = F,
       xlab = "", ylab = "", ...)

  graphics::title(xlab="Estimated label", line=0)
  graphics::title(ylab="True label", line=0)

  graphics::rect(0, 0, 1, 1)
  graphics::rect(0, .5, .5 , 1, col = grDevices::rgb(col1[1], col1[2], col1[3]))
  graphics::rect(.5, 0, 1 , .5, col = grDevices::rgb(col2[1], col2[2], col2[3]))

  graphics::text(0.25,0.75, labels = paste0("Dropout,\nn = ", tab[1,1], "\n(", round(prob1*100), "%)"))
  graphics::text(0.75,0.75, labels = paste0("n = ", tab[1,2], "\n(", 100-round(prob1*100), "%)"))
  graphics::text(0.75,0.25, labels = paste0("Not dropout, \nn = ", tab[2,2], "\n(", round(prob2*100), "%)"))
  graphics::text(0.25,0.25, labels = paste0("n = ", tab[2,1], "\n(", 100-round(prob2*100), "%)"))

  invisible()
}
