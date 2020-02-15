rm(list=ls())

load("../results/factorization_esvd.RData")

# clean results
for(i in 1:length(res)){
  res[[i]] <- res[[i]][which(sapply(res[[i]], function(x){!all(is.na(x))}))]
}
sapply(res, function(x){
  length(x)
})


# plot the embeddings
k <- 1
par(mfrow = c(3,3))
for(i in 1:8){
  plot(res[[(k-1)*8+i]][[1]]$fit[,1], res[[(k-1)*8+i]][[1]]$fit[,2], asp = T, col = rep(1:4, each = paramMat[1,"n_each"]),
       pch = 16, main = i)
}

# first plot the "true values" and its observed values, as a sanity check
correct_idx <- c(1, 10, 20, 31)
par(mfrow = c(2,2))
for(i in correct_idx){
  plot_mat <- lapply(1:length(res[[i]]), function(j){
    cbind(res[[i]][[j]]$missing_val, res[[i]][[j]]$expected_val)
  })

  if(length(plot_mat) > 1) plot_mat <- do.call(rbind, plot_mat) else plot_mat <- plot_mat[[1]]

  plot(plot_mat[,2], plot_mat[,1], asp = T, main = i, pch = 16, xlab = "True expected value", ylab = "Observed value")
  lines(c(-1e5,1e5), c(-1e5, 1e5), col = "red", lty = 2, lwd = 2)

  seq_val <- seq(0, 4000, length.out = 500)
  if(i == 1){
    # estimate std
    sd_val <- stats::sd(plot_mat[,1] - plot_mat[,2])

    y_bot <- sapply(seq_val, function(x){
      stats::qnorm(0.1, mean = x, sd = sd_val)
    })
    y_top <- sapply(seq_val, function(x){
      stats::qnorm(0.9, mean = x, sd = sd_val)
    })
  } else if(i == 10){
    y_bot <- sapply(seq_val, function(x){
      stats::qpois(0.1, lambda = x)
    })
    y_top <- sapply(seq_val, function(x){
      stats::qpois(0.9, lambda = x)
    })
  } else if(i == 20){
    size <- as.numeric(paramMat[i,"fitting_param"])
    y_bot <- sapply(seq_val, function(x){
      stats::qnbinom(0.1, size = size, prob = size/(size+x))
    })
    y_top <- sapply(seq_val, function(x){
      stats::qnbinom(0.9, size = size, prob = size/(size+x))
    })
#
#     y_bot <- sapply(seq_val, function(x){
#       p_val <- 1-size/(size+x)
#       x - 0.67*(size*p_val)/(p_val)^2
#     })
#     y_top <- sapply(seq_val, function(x){
#       p_val <- 1-size/(size+x)
#       x + 0.67*(size*p_val)/(p_val)^2
#     })
#     y_bot[is.nan(y_bot)] <- 0
#     y_top[is.nan(y_top)] <- 0

  } else {
    scalar <- as.numeric(paramMat[i,"fitting_param"])
    y_bot <- sapply(seq_val, function(x){
      stats::qnorm(0.1, mean = x, sd = x/scalar)
    })
    y_top <- sapply(seq_val, function(x){
      stats::qnorm(0.9, mean = x, sd = x/scalar)
    })
  }

  lines(seq_val, y_bot, col = "red", lty = 1, lwd = 2)
  lines(seq_val, y_top, col = "red", lty = 1, lwd = 2)

  pca_res <- stats::princomp(plot_mat)
  lines(c(-1e10,1e10)*pca_res$loadings[1,2], c(-1e10,1e10)*pca_res$loadings[1,1],
        col = "blue", lwd = 2, lty = 2)
}

###############################
#
# # plot the missing values
correct_idx <- c(1, 10, 20, 31)
par(mfrow = c(2,2))
for(i in correct_idx){
  plot_mat <- lapply(1:length(res[[i]]), function(j){
    cbind(res[[i]][[j]]$missing_val, res[[i]][[j]]$pred_val)
  })

  if(length(plot_mat) > 1) plot_mat <- do.call(rbind, plot_mat) else plot_mat <- plot_mat[[1]]

  plot(plot_mat[,2], plot_mat[,1], asp = T, main = i, pch = 16, xlab = "True expected value", ylab = "Observed value")
  lines(c(-1e5,1e5), c(-1e5, 1e5), col = "red", lty = 2, lwd = 2)

  seq_val <- seq(0, 4000, length.out = 500)
  if(i == 1){
    # estimate std
    sd_val <- stats::sd(plot_mat[,1] - plot_mat[,2])

    y_bot <- sapply(seq_val, function(x){
      stats::qnorm(0.1, mean = x, sd = sd_val)
    })
    y_top <- sapply(seq_val, function(x){
      stats::qnorm(0.9, mean = x, sd = sd_val)
    })
  } else if(i == 10){
    y_bot <- sapply(seq_val, function(x){
      stats::qpois(0.1, lambda = x)
    })
    y_top <- sapply(seq_val, function(x){
      stats::qpois(0.9, lambda = x)
    })
  } else if(i == 20){
    size <- as.numeric(paramMat[i,"fitting_param"])
    y_bot <- sapply(seq_val, function(x){
      stats::qnbinom(0.1, size = size, prob = size/(size+x))
    })
    y_top <- sapply(seq_val, function(x){
      stats::qnbinom(0.9, size = size, prob = size/(size+x))
    })
  } else {
    scalar <- as.numeric(paramMat[i,"fitting_param"])
    y_bot <- sapply(seq_val, function(x){
      stats::qnorm(0.1, mean = x, sd = x/scalar)
    })
    y_top <- sapply(seq_val, function(x){
      stats::qnorm(0.9, mean = x, sd = x/scalar)
    })
  }

  lines(seq_val, y_bot, col = "red", lty = 1, lwd = 2)
  lines(seq_val, y_top, col = "red", lty = 1, lwd = 2)

  pca_res <- stats::princomp(plot_mat)
  lines(c(-1e10,1e10)*pca_res$loadings[1,2], c(-1e10,1e10)*pca_res$loadings[1,1],
        col = "blue", lwd = 2, lty = 2)
}




# k <- 1
# par(mfrow = c(3,3))
# for(i in 1:8){
#   plot_mat <- lapply(1:length(res[[(k-1)*8+i]]), function(j){
#     cbind(res[[(k-1)*8+i]][[j]]$missing_val, res[[(k-1)*8+i]][[j]]$pred_val)
#   })
#
#   if(length(plot_mat) > 1) plot_mat <- do.call(rbind, plot_mat) else plot_mat <- plot_mat[[1]]
#
#   plot(plot_mat[,2], plot_mat[,1], asp = T, main = i, pch = 16)
#   lines(c(-1e5,1e5), c(-1e5, 1e5), col = "red", lty = 2, lwd = 2)
#
#   seq_val <- seq(0, 4000, length.out = 500)
#   if(i == 1){
#     # estimate std
#     sd_val <- stats::sd(plot_mat[,1] - plot_mat[,2])
#
#     y_bot <- sapply(seq_val, function(x){
#       stats::qnorm(0.1, mean = x, sd = sd_val)
#     })
#     y_top <- sapply(seq_val, function(x){
#       stats::qnorm(0.9, mean = x, sd = sd_val)
#     })
#   } else if(i == 2){
#     y_bot <- sapply(seq_val, function(x){
#       stats::qpois(0.1, lambda = x)
#     })
#     y_top <- sapply(seq_val, function(x){
#       stats::qpois(0.9, lambda = x)
#     })
#   } else if(i %in% c(3:5)){
#     size <- as.numeric(paramMat[i,"fitting_param"])
#     # y_bot <- sapply(seq_val, function(x){
#     #   stats::qnbinom(0.25, size = size, prob = 1-size/(size+x))
#     # })
#     # y_top <- sapply(seq_val, function(x){
#     #   stats::qnbinom(0.75, size = size, prob = 1-size/(size+x))
#     # })
#
#     y_bot <- sapply(seq_val, function(x){
#       p_val <- 1-size/(size+x)
#       x - 0.67*(size*p_val)/(p_val)^2
#     })
#     y_top <- sapply(seq_val, function(x){
#       p_val <- 1-size/(size+x)
#       x + 0.67*(size*p_val)/(p_val)^2
#     })
#     y_bot[is.nan(y_bot)] <- 0
#     y_top[is.nan(y_top)] <- 0
#
#   } else {
#     scalar <- as.numeric(paramMat[i,"fitting_param"])
#     y_bot <- sapply(seq_val, function(x){
#       stats::qnorm(0.25, mean = x, sd = x/scalar)
#     })
#     y_top <- sapply(seq_val, function(x){
#       stats::qnorm(0.75, mean = x, sd = x/scalar)
#     })
#   }
#
#   lines(seq_val, y_bot, col = "red", lty = 1, lwd = 2)
#   lines(seq_val, y_top, col = "red", lty = 1, lwd = 2)
#
#   pca_res <- stats::princomp(plot_mat)
#   lines(c(-1e10,1e10)*pca_res$loadings[1,2], c(-1e10,1e10)*pca_res$loadings[1,1],
#         col = "blue", lwd = 2, lty = 2)
# }
