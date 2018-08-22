rm(list=ls())
source("../experiment/em_gamma_normal.R")
source("../experiment/em_gamma_truncatednormal.R")
source("../experiment/truncated_normal_estimation.R")
source("../experiment/truncated_normal_estimation_MoM.R")

generate_univariate_data <- function(n, vec, min_val = log10(1.01),
                                     return_assignment = F){
  prop = vec[1]; shape = vec[2]; rate = vec[3]; mu = vec[4]; sd = vec[5]

  assign_vec <- sample(2, n, replace = T, prob = c(prop, 1-prop))
  x <- numeric(n)

  idx1 <- which(assign_vec == 1)
  x[idx1] <- stats::rgamma(length(idx1), shape = shape, rate = rate)

  idx2 <- which(assign_vec == 2)
  x[idx2] <- stats::rnorm(length(idx2), mean = mu, sd = sd)

  while(sum(x < min_val) > 0){
    idx_tmp <- which(x < min_val)
    if(length(idx_tmp) == 0) break()
    x[idx_tmp] <- stats::rnorm(length(idx_tmp), mean = mu, sd = sd)
  }

  if(return_assignment){
    list(x = x, assign_vec = assign_vec)
  } else {
    x
  }
}


draw_curve <- function(vec_list, ymax = NA, max_val = 15, min_val = log10(1.01)/10, ...){
  #compute distribution
  x_seq <- seq(min_val, max_val, length.out = 1000)

  density_list <- lapply(vec_list, function(vec){
    gamma_dist <- vec[1] * sapply(x_seq, stats::dgamma, shape = vec[2], rate = vec[3])
    norm_dist <- (1-vec[1]) * sapply(x_seq, stats::dnorm, mean = vec[4], sd = vec[5])
    norm_dist <- norm_dist / (1 - stats::pnorm(0, mean = vec[4], sd = vec[5]))

    list(gamma_dist, norm_dist)
  })

  if(is.na(ymax)) ymax <- max(unlist(density_list))

  plot(NA, xlim = c(min_val, max_val), ylim = c(0, ymax), ...)
  for(i in 1:length(vec_list)){
    lines(x_seq, density_list[[i]][[1]], col = rgb(0.584, 0.858, 0.564), lty = i, ... )
    lines(x_seq, density_list[[i]][[2]], col = rgb(0.803, 0.156, 0.211), lty = i, ...)
  }

  invisible()
}

compute_error <- function(vec1, vec2, interval = 0.01){
  x_max <- max(vec1[4]+5*vec1[5], vec2[4]+5*vec2[5])
  x_seq <- seq(0, x_max, by = interval)

  norm_dist1 <- (1-vec1[1]) * sapply(x_seq, stats::dnorm, mean = vec1[4], sd = vec1[5])
  norm_dist1 <- norm_dist1 / (1 - stats::pnorm(0, mean = vec1[4], sd = vec1[5]))

  norm_dist2 <- (1-vec2[1]) * sapply(x_seq, stats::dnorm, mean = vec2[4], sd = vec2[5])
  norm_dist2 <- norm_dist2 / (1 - stats::pnorm(0, mean = vec2[4], sd = vec2[5]))

  res <- c(sum(interval * (norm_dist1 - norm_dist2)^2), max(abs(norm_dist1 - norm_dist2)),
           abs(vec1[4] - vec2[4]))
  names(res) <- c("l2", "l0", "mean")

  res
}


######################

set.seed(10)
res_pop <- c(.8, 1, 10, 1.5, 1)
x <- generate_univariate_data(1000, res_pop)
res_trunc_8 <- .get_mix(x, MoM = T, prop_init = .8)
res_trunc_9 <- .get_mix(x, MoM = T, prop_init = .9)

err_8 <- compute_error(res_pop, res_trunc_8)
err_9 <- compute_error(res_pop, res_trunc_9)

png("../figure/experiment/17_experiment_mom.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))
draw_curve(list(res_pop, res_trunc_8), lwd = 3, max_val = 3, ymax = 0.5,
           main = paste0("MoM, Initialization 1:\nError = (", paste0(round(err_8, 2),
                                                       collapse = ", "), ")"),
           xlab = "Value", ylab = "Density")
draw_curve(list(res_pop, res_trunc_9), lwd = 3, max_val = 3, ymax = 0.5,
           main = paste0("MoM, Initialization 2:\nError = (", paste0(round(err_9, 2),
                                                        collapse = ", "), ")"),
           xlab = "Value", ylab = "Density")
graphics.off()

