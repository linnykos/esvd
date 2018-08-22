estimate_truncated_normal_mom <- function(vec, weight = rep(1, length(vec)),
                                      breaks = 100, iterations = 5, window = 0.3){
  mean_val <- sum(weight * vec)/sum(weight)
  var_val <- sum(weight * (vec - mean_val)^2)/sum(weight)

  m_range <- c(0, max(vec))
  v_range <- c(0.1, 2*diff(range(vec)))
  m_est <- mean_val; v_est <- var_val
  window_size <- ceiling(window * breaks)

  for(i in 1:iterations){
    m_vec <- seq(m_range[1], m_range[2], length.out = breaks)
    diff_vec <- compute_difference_vector(m_vec, breaks, mean_val,
                                          var_val, sigma = sqrt(v_est))
    idx <- which.min(diff_vec)
    m_est <- m_vec[idx]
    m_range <- m_vec[c(max(idx-window_size, 1), min(idx+window_size, breaks))]

    v_vec <- seq(v_range[1], v_range[2], length.out = breaks)
    diff_vec <- compute_difference_vector(v_vec, breaks, mean_val,
                                          var_val, mu = m_est)
    idx <- which.min(diff_vec)
    v_est <- v_vec[idx]
    v_range <- v_vec[c(max(idx-window_size, 1), min(idx+window_size, breaks))]
  }

  list(mean = m_est, var = v_est)
}

compute_difference_vector <- function(val_vec, breaks, mean_val, var_val,
                                      ...){
  mean_tmp <- sapply(val_vec, compute_mean, ...)
  var_tmp <- sapply(val_vec, compute_var, ...)

  compute_difference(mean_tmp, var_tmp, mean_val, var_val)
}

.l2norm <- function(x){sqrt(sum(x^2))}

compute_difference <- function(vec1, vec2, target1, target2){
  mat <- cbind(vec1 - target1, vec2 - target2)
  apply(mat, 1, .l2norm)
}

compute_mean <- function(mu, sigma){
  mu + sigma * stats::dnorm(0)/(1-stats::pnorm(0))
}

compute_var <- function(mu, sigma){
  alpha <- (0-mu)/sigma
  Z <- 1-stats::pnorm(0)
  sigma^2 * (1+ alpha*stats::dnorm(alpha)/Z - (stats::dnorm(alpha)/Z)^2)
}
