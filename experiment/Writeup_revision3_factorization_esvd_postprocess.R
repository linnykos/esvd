rm(list=ls())

load("../results/factorization_esvd.RData")

# clean results
for(i in 1:length(res)){
  res[[i]] <- res[[i]][which(sapply(res[[i]], function(x){!all(is.na(x))}))]
}
table(as.numeric(sapply(res, function(x){
  length(x)
})))

################################

correct_idx <- c(1, 10, 20, 31)
main_vec <- c("Gaussian\n(constant variance)", "Poisson",
              "Negative binomial\n(size = 50)", "Curved Gaussian\n(alpha = 2)")
main_vec2 <- c("Gaussian (constant variance)", "Poisson",
               "Negative binomial (size = 50)", "Curved Gaussian (alpha = 2)")

png(filename = "../figure/experiment/Revision_writeup3_simulation_missing_value.png",
    height = 2250, width = 2250, res = 300,
    units = "px")
par(mfrow = c(2,2))
for(k in 1:4){
  i <- correct_idx[k]
  plot_mat <- lapply(1:length(res[[i]]), function(j){
    cbind(res[[i]][[j]]$missing_val, res[[i]][[j]]$pred_val)
  })

  if(length(plot_mat) > 1) plot_mat <- do.call(rbind, plot_mat) else plot_mat <- plot_mat[[1]]

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
    size <- as.numeric(paramMat[i,"true_r"])
    y_bot <- sapply(seq_val, function(x){
      stats::qnbinom(0.1, size = size, prob = size/(size+x))
    })
    y_top <- sapply(seq_val, function(x){
      stats::qnbinom(0.9, size = size, prob = size/(size+x))
    })
  } else {
    scalar <- as.numeric(paramMat[i,"true_scalar"])
    y_bot <- sapply(seq_val, function(x){
      stats::qnorm(0.1, mean = x, sd = x/scalar)
    })
    y_top <- sapply(seq_val, function(x){
      stats::qnorm(0.9, mean = x, sd = x/scalar)
    })
  }


  plot(NA, asp = T, main = main_vec[k], xlim = range(plot_mat[,2]), ylim = range(plot_mat[,1]),
       xlab = "Predicted missing value", ylab = "Observed value")

  polygon(c(seq_val, rev(seq_val)), c(y_top, rev(y_bot)), col = rgb(1,0,0,0.2),
          border = NA, density = 30, angle = -45)
  points(plot_mat[,2], plot_mat[,1], pch = 16, col = rgb(0,0,0,0.2))

  lines(rep(0,2), c(-1e10,1e10), col = "red", lwd = 1)
  lines(c(-1e10,1e10), rep(0,2), col = "red", lwd = 1)
  lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)

  lines(seq_val, y_bot, col = "red", lty = 2, lwd = 2)
  lines(seq_val, y_top, col = "red", lty = 2, lwd = 2)

  pca_res <- stats::prcomp(plot_mat, center = F, scale = F)
  lines(c(0, 1e6), c(0, 1e6*pca_res$rotation[2,1]/pca_res$rotation[1,1]), col = "blue", lwd = 2, lty = 2)
  #
  #   lines(c(-1e10,1e10)*pca_res$rotation[1,2], c(-1e10,1e10)*pca_res$rotation[1,1],
  #         col = "blue", lwd = 2, lty = 2)

  rad <- 2/5*max(plot_mat[,1])
  ang <- as.numeric(acos(abs(c(1,0) %*% pca_res$rotation[,1])))
  radian_seq <- seq(0, ang, length.out = 100)
  x_circ <- rad * cos(radian_seq)
  y_circ <- rad * sin(radian_seq)
  lines(x_circ, y_circ, lty = 2)
  text(x = rad , y = 2/5*rad, pos = 4, label = paste0(round(ang * 180/pi, 1), " degrees"))
}

graphics.off()

####################################

png(filename = "../figure/experiment/Revision_writeup3_simulation_fit.png",
    height = 2250, width = 2250, res = 300,
    units = "px")
pos_vec <- c("topleft", "topleft", "topright", "topright")
par(mfrow = c(2,2))
for(k in 1:4){
  i <- correct_idx[k]
  plot(res[[i]][[1]]$fit[,1], res[[i]][[1]]$fit[,2], col = rep(1:4, each = paramMat[1,"n_each"]), pch = 16,
       xlab = "Latent dimension 1", ylab = "Latent dimension 2", main =
         paste0("Estimated cell latent positions\n(", main_vec2[k], ")"), asp = T)
  legend(pos_vec[k], c("Cell type 1", "Cell type 2", "Cell type 3", "Cell type 4"),
         bty = "n", fill = 1:4)
}
graphics.off()

################################

correct_idx <- c(1, 2, 4, 7)
main_vec <- c("Gaussian (fixed variance)", "Poisson",
              "Neg. bin. (size = 25)",
              "Neg. bin. (size = 50)",
              "Neg. bin. (size = 200)",
              "Curved Gaussian (scalar = 1)",
              "Curved Gaussian (scalar = 2)",
              "Curved Gaussian (scalar = 4)",
              "Exponential")

for(k in 1:4){
  png(filename = paste0("../figure/experiment/Revision_writeup3_simulation_missing_", k, ".png"),
      height = 2500, width = 2500, res = 300,
      units = "px")

  par(mfrow = c(3,3), mar = c(4,4,4,0.5))
  for(i in 1:9){
    plot_mat <- lapply(1:length(res[[(k-1)*8+i]]), function(j){
      cbind(res[[(k-1)*8+i]][[j]]$missing_val, res[[(k-1)*8+i]][[j]]$pred_val)
    })

    if(length(plot_mat) > 1) plot_mat <- do.call(rbind, plot_mat) else plot_mat <- plot_mat[[1]]

    seq_val <- seq(0, 4000, length.out = 500)
    if(i == 1){
      sd_val <- stats::sd(plot_mat[,1] - plot_mat[,2])

      y_bot <- sapply(seq_val, function(x){
        stats::qnorm(0.1, mean = x, sd = sd_val)
      })
      y_top <- sapply(seq_val, function(x){
        stats::qnorm(0.9, mean = x, sd = sd_val)
      })
    } else if(i == 2){
      y_bot <- sapply(seq_val, function(x){
        stats::qpois(0.1, lambda = x)
      })
      y_top <- sapply(seq_val, function(x){
        stats::qpois(0.9, lambda = x)
      })
    } else if(i %in% c(3:5)){
      size <- as.numeric(paramMat[i,"fitting_param"])
      y_bot <- sapply(seq_val, function(x){
        stats::qnbinom(0.1, size = size, prob = size/(size+x))
      })
      y_top <- sapply(seq_val, function(x){
        stats::qnbinom(0.9, size = size, prob = size/(size+x))
      })
    } else if(i %in% (6:8)){
      scalar <- as.numeric(paramMat[i,"fitting_param"])
      y_bot <- sapply(seq_val, function(x){
        stats::qnorm(0.1, mean = x, sd = x/scalar)
      })
      y_top <- sapply(seq_val, function(x){
        stats::qnorm(0.9, mean = x, sd = x/scalar)
      })
    } else {
      y_bot <- sapply(seq_val, function(x){
        stats::qexp(0.1, rate = 1/x)
      })
      y_top <- sapply(seq_val, function(x){
        stats::qexp(0.9, rate = 1/x)
      })
    }

    plot(NA, asp = T, main = paste0(main_vec[i], ifelse(i == correct_idx[k], "\n(True model)", "")), pch = 16,
         xlim = range(plot_mat[,2]), ylim = range(plot_mat[,1]),
         xlab = "Predicted missing value", ylab = "Observed value")

    polygon(c(seq_val, rev(seq_val)), c(y_top, rev(y_bot)), col = rgb(1,0,0,0.2),
            border = NA, density = 30, angle = -45)
    points(plot_mat[,2], plot_mat[,1], pch = 16, col = rgb(0,0,0,0.2))

    lines(rep(0,2), c(-1e10,1e10), col = "red", lwd = 1)
    lines(c(-1e10,1e10), rep(0,2), col = "red", lwd = 1)
    lines(c(-1e10,1e10), c(-1e10,1e10), col = "red", lwd = 2)

    lines(seq_val, y_bot, col = "red", lty = 2, lwd = 2)
    lines(seq_val, y_top, col = "red", lty = 2, lwd = 2)

    pca_res <- stats::prcomp(plot_mat, center = F, scale = F)
    lines(c(0, 1e6), c(0, 1e6*pca_res$rotation[2,1]/pca_res$rotation[1,1]), col = "blue", lwd = 2, lty = 2)

  }

  graphics.off()
}

#################################

correct_idx <- c(1, 2, 4, 7)
main_vec <- c("Gaussian (fixed variance)", "Poisson",
              "Neg. bin. (size = 25)",
              "Neg. bin. (size = 50)",
              "Neg. bin. (size = 200)",
              "Curved Gaussian (scalar = 1)",
              "Curved Gaussian (scalar = 2)",
              "Curved Gaussian (scalar = 4)",
              "Exponential")

for(k in 1:4){
  png(filename = paste0("../figure/experiment/Revision_writeup3_simulation_embedding_", k, ".png"),
      height = 2500, width = 2500, res = 300,
      units = "px")

  par(mfrow = c(3,3), mar = c(4,4,4,0.5))
  for(i in 1:9){
    plot(res[[(k-1)*8+i]][[1]]$fit[,1], res[[(k-1)*8+i]][[1]]$fit[,2], col = rep(1:4, each = paramMat[1,"n_each"]), pch = 16,
         xlab = "Latent dimension 1", ylab = "Latent dimension 2", main =
           paste0(main_vec[i], ifelse(i == correct_idx[k], "\n(True model)", "")), asp = T)
  }

  graphics.off()
}

