estimate_truncated_normal <- function(vec, breaks = 100, iterations = 5, window = 0.3){
  mean_val <- mean(vec); var_val <- var(vec)

  m_range <- c(-1,1)*max(abs(vec))
  v_range <- c(0.1, 2*diff(range(vec)))
  m_est <- mean_val; v_est <- var_val
  window_size <- ceiling(window * breaks)

  for(i in 1:iterations){
    diff_vec <- compute_difference_vector(vec, m_range, breaks, mean_val,
                                          var_val, sigma = sqrt(v_est))
    idx <- which.min(diff_vec)
    m_est <- m_idx[idx]
    m_range <- m_idx[c(max(idx-window_size, 1), min(idx+window_size, breaks))]

    diff_vec <- compute_difference_vector(vec, m_range, breaks, mean_val,
                                          var_val, mu = m_est)
    idx <- which.min(diff_vec)
    v_est <- v_idx[idx]
    v_range <- v_idx[c(max(idx-window_size, 1), min(idx+window_size, breaks))]
  }

  list(mean = m_est, var = v_est)
}

compute_difference_vector <- function(vec, range_vec, breaks, mean_val, var_val,
                                      ...){
  val_vec <- seq(range_vec[1], range_vec[2], length.out = breaks)
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
