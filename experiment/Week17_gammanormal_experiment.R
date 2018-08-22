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

  res <- c(sum(interval * (norm_dist1 - norm_dist2)^2), max(norm_dist1 - norm_dist2),
           abs(vec1[4] - vec2[4]))
  names(res) <- c("l2", "l0", "mean")

  res
}

compute_dropout_table <- function(data, vec2){
  gamma_dist <- vec2[1] * sapply(data$x, stats::dgamma, shape = vec2[2], rate = vec2[3])
  norm_dist <- (1-vec2[1]) * sapply(data$x, stats::dnorm, mean = vec2[4], sd = vec2[5])
  norm_dist <- norm_dist / (1 - stats::pnorm(0, mean = vec2[4], sd = vec2[5]))

  ratio <- gamma_dist/(gamma_dist + norm_dist)
  estimated <- rep(2, length(ratio))
  estimated[ratio > 0.5] <- 1

  truth <- data$assign_vec

  table(truth, estimated)
}

visualize_table <- function(data, vec, ...){
  tab <- compute_dropout_table(data, vec)

  col1 <- c(0.584, 0.858, 0.564)
  col2 <- c(0.803, 0.156, 0.211)

  prob1 <- tab[1,1]/sum(tab[1,])
  prob2 <- tab[2,2]/sum(tab[2,])

  col1 <- sapply(col1, function(x){1 - prob1 * (1 - x)})
  col2 <- sapply(col2, function(x){1 - prob2 * (1 - x)})

  par(mar = c(1,1,3,0))
  plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T, axes = F,
       xlab = "", ylab = "", ...)

  title(xlab="Estimated label", line=0)
  title(ylab="True label", line=0)

  rect(0, 0, 1, 1)
  rect(0, .5, .5 , 1, col = rgb(col1[1], col1[2], col1[3]))
  rect(.5, 0, 1 , .5, col = rgb(col2[1], col2[2], col2[3]))

  text(0.25,0.75, labels = paste0("Dropout,\nn = ", tab[1,1], "\n(", round(prob1*100), "%)"))
  text(0.75,0.75, labels = paste0("n = ", tab[1,2], "\n(", 100-round(prob1*100), "%)"))
  text(0.75,0.25, labels = paste0("Not dropout, \nn = ", tab[2,2], "\n(", round(prob2*100), "%)"))
  text(0.25,0.25, labels = paste0("n = ", tab[2,1], "\n(", 100-round(prob2*100), "%)"))

  invisible()
}



###

set.seed(10)
res_pop <- c(.1, .5, 10, .5, 1)
x <- generate_univariate_data(1000, res_pop)
res_original <- get_mix(x, prop_init = .3)
res_trunc <- .get_mix(x, prop_init = .3)
# res_trunc_MoM <- .get_mix(x, MoM = T, prop_init = .3)

err_original <- compute_error(res_pop, res_original)
err_trunc <- compute_error(res_pop, res_trunc)

png("../figure/experiment/17_experiment_density_small_near.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))
draw_curve(list(res_pop, res_original), lwd = 3, max_val = 3, ymax = 2,
           main = paste0("Original: Error = (", paste0(round(err_original, 2),
                                               collapse = ", "), ")"),
           xlab = "Value", ylab = "Density")
draw_curve(list(res_pop, res_trunc), lwd = 3, max_val = 3, ymax = 2,
           main = paste0("Truncated: Error = (", paste0(round(err_trunc, 2),
                                                       collapse = ", "), ")"),
           xlab = "Value", ylab = "Density")
graphics.off()

set.seed(10)
new_x <- generate_univariate_data(1000, res_pop, return_assignment = T)
png("../figure/experiment/17_experiment_density_small_near_table.png", height = 800, width = 1800, res = 300, units = "px")
par(mfrow = c(1,3))
visualize_table(new_x, res_pop, main = "True parameters")
visualize_table(new_x, res_original, main = "Estimated parameters\n(original)")
visualize_table(new_x, res_trunc, main = "Estimated parameters\n(truncated)")
graphics.off()

###

set.seed(10)
res_pop <- c(.3, .5, 10, .5, 1)
x <- generate_univariate_data(1000, res_pop)
res_original <- get_mix(x, prop_init = .3)
res_trunc <- .get_mix(x, prop_init = .3)
# res_trunc_MoM <- .get_mix(x, MoM = T, prop_init = .3)

err_original <- compute_error(res_pop, res_original)
err_trunc <- compute_error(res_pop, res_trunc)

png("../figure/experiment/17_experiment_density_medium_near.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))
draw_curve(list(res_pop, res_original), lwd = 3, max_val = 3, ymax = 2,
           main = paste0("Original: Error = (", paste0(round(err_original, 2),
                                                       collapse = ", "), ")"),
           xlab = "Value", ylab = "Density")
draw_curve(list(res_pop, res_trunc), lwd = 3, max_val = 3, ymax = 2,
           main = paste0("Truncated: Error = (", paste0(round(err_trunc, 2),
                                                        collapse = ", "), ")"),
           xlab = "Value", ylab = "Density")
graphics.off()

set.seed(10)
new_x <- generate_univariate_data(1000, res_pop, return_assignment = T)
png("../figure/experiment/17_experiment_density_medium_near_table.png", height = 800, width = 1800, res = 300, units = "px")
par(mfrow = c(1,3))
visualize_table(new_x, res_pop, main = "True parameters")
visualize_table(new_x, res_original, main = "Estimated parameters\n(original)")
visualize_table(new_x, res_trunc, main = "Estimated parameters\n(truncated)")
graphics.off()

###

set.seed(10)
res_pop <- c(.3, 1, 1, 1.5, 1)
x <- generate_univariate_data(1000, res_pop)
res_original <- get_mix(x, prop_init = .3)
res_trunc <- .get_mix(x, prop_init = .3)
# res_trunc_MoM <- .get_mix(x, MoM = T, prop_init = .3)

err_original <- compute_error(res_pop, res_original)
err_trunc <- compute_error(res_pop, res_trunc)

png("../figure/experiment/17_experiment_density_medium_far.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))
draw_curve(list(res_pop, res_original), lwd = 3, max_val = 3, ymax = 0.5,
           main = paste0("Original: Error = (", paste0(round(err_original, 2),
                                                       collapse = ", "), ")"),
           xlab = "Value", ylab = "Density")
draw_curve(list(res_pop, res_trunc), lwd = 3, max_val = 3, ymax = 0.5,
           main = paste0("Truncated: Error = (", paste0(round(err_trunc, 2),
                                                        collapse = ", "), ")"),
           xlab = "Value", ylab = "Density")
graphics.off()

set.seed(10)
new_x <- generate_univariate_data(1000, res_pop, return_assignment = T)
png("../figure/experiment/17_experiment_density_medium_far_table.png", height = 800, width = 1800, res = 300, units = "px")
par(mfrow = c(1,3))
visualize_table(new_x, res_pop, main = "True parameters")
visualize_table(new_x, res_original, main = "Estimated parameters\n(original)")
visualize_table(new_x, res_trunc, main = "Estimated parameters\n(truncated)")
graphics.off()

###

set.seed(10)
res_pop <- c(.1, 1, 1, 1.5, 1)
x <- generate_univariate_data(1000, res_pop)
res_original <- get_mix(x, prop_init = .3)
res_trunc <- .get_mix(x, prop_init = .3)
# res_trunc_MoM <- .get_mix(x, MoM = T, prop_init = .3)

err_original <- compute_error(res_pop, res_original)
err_trunc <- compute_error(res_pop, res_trunc)

png("../figure/experiment/17_experiment_density_small_far.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))
draw_curve(list(res_pop, res_original), lwd = 3, max_val = 3, ymax = 2,
           main = paste0("Original: Error = (", paste0(round(err_original, 2),
                                                       collapse = ", "), ")"),
           xlab = "Value", ylab = "Density")
draw_curve(list(res_pop, res_trunc), lwd = 3, max_val = 3, ymax = 2,
           main = paste0("Truncated: Error = (", paste0(round(err_trunc, 2),
                                                        collapse = ", "), ")"),
           xlab = "Value", ylab = "Density")
graphics.off()

set.seed(10)
new_x <- generate_univariate_data(1000, res_pop, return_assignment = T)
png("../figure/experiment/17_experiment_density_small_far_table.png", height = 800, width = 1800, res = 300, units = "px")
par(mfrow = c(1,3))
visualize_table(new_x, res_pop, main = "True parameters")
visualize_table(new_x, res_original, main = "Estimated parameters\n(original)")
visualize_table(new_x, res_trunc, main = "Estimated parameters\n(truncated)")
graphics.off()

###

set.seed(10)
res_pop <- c(.3, 1, 1, 4, 1)
x <- generate_univariate_data(1000, res_pop)
res_original <- get_mix(x, prop_init = .3)
res_trunc <- .get_mix(x, prop_init = .3)
# res_trunc_MoM <- .get_mix(x, MoM = T, prop_init = .3)

err_original <- compute_error(res_pop, res_original)
err_trunc <- compute_error(res_pop, res_trunc)

png("../figure/experiment/17_experiment_density_medium_reallyfar.png", height = 1200, width = 2400, res = 300, units = "px")
par(mfrow = c(1,2))
draw_curve(list(res_pop, res_original), lwd = 3, ymax = .5, max_val = 7,
           main = paste0("Original: Error = (", paste0(round(err_original, 2),
                                                       collapse = ", "), ")"),
           xlab = "Value", ylab = "Density")
draw_curve(list(res_pop, res_trunc), lwd = 3, ymax = .5, max_val = 7,
           main = paste0("Truncated: Error = (", paste0(round(err_trunc, 2),
                                                        collapse = ", "), ")"),
           xlab = "Value", ylab = "Density")
graphics.off()

set.seed(10)
new_x <- generate_univariate_data(1000, res_pop, return_assignment = T)
png("../figure/experiment/17_experiment_density_medium_reallyfar_table.png", height = 800, width = 1800, res = 300, units = "px")
par(mfrow = c(1,3))
visualize_table(new_x, res_pop, main = "True parameters")
visualize_table(new_x, res_original, main = "Estimated parameters\n(original)")
visualize_table(new_x, res_trunc, main = "Estimated parameters\n(truncated)")
graphics.off()

