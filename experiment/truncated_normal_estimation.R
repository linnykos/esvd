estimate_truncated_normal <- function(vec, weight = rep(1, length(vec)),
                                      min_val = log10(1.01), multiplier = 5,
                                      breaks = 100, iterations = 5, window = 0.3){
  n <- length(vec)
  mean_val <- sum(weight * vec)/sum(weight)
  var_val <- sum(weight * (vec - mean_val)^2)/sum(weight)

  obs_lis <- .estimate_cdf(vec, weight)

  m_range <- c(0, max(vec))
  v_range <- c(0.1, 2*diff(range(vec)))
  m_est <- mean_val; v_est <- var_val
  window_size <- ceiling(window * breaks)

  for(i in 1:iterations){
    m_vec <- seq(m_range[1], m_range[2], length.out = breaks)
    diff_vec <- sapply(m_vec, function(x){
      tmp <- stats::rnorm(multiplier*n, mean = x, sd = sqrt(v_est))
      tmp <- tmp[tmp > min_val]
      tmp_lis <- .estimate_cdf(tmp)
      .ks_distance(obs_lis, tmp_lis)
    })

    idx <- which.min(diff_vec)
    m_est <- m_vec[idx]
    m_range <- m_vec[c(max(idx-window_size, 1), min(idx+window_size, breaks))]

    v_vec <- seq(v_range[1], v_range[2], length.out = breaks)
    diff_vec <- sapply(v_vec, function(x){
      tmp <- stats::rnorm(multiplier*n, mean = m_est, sd = x)
      tmp <- tmp[tmp > min_val]
      tmp_lis <- .estimate_cdf(tmp)
      .ks_distance(obs_lis, tmp_lis)
    })

    idx <- which.min(diff_vec)
    v_est <- v_vec[idx]
    v_range <- v_vec[c(max(idx-window_size, 1), min(idx+window_size, breaks))]
  }

  list(mean = m_est, var = v_est)
}

.estimate_cdf <- function(vec, weight = rep(1, length(vec))){
  idx <- order(vec)
  vec <- vec[idx]; weight <- weight[idx]

  uniq_val <- sort(unique(vec))
  cdf_vec <- cumsum(weight)/sum(weight)

  list(x = uniq_val, y = cdf_vec)
}

.ks_distance <- function(lis1, lis2){
  all_x <- c(lis1$x, lis2$x)
  all_y <- c(lis1$y, lis2$y)
  idx <- order(all_x)
  all_x <- all_x[idx]; all_y <- all_y[idx]
  n <- length(all_y)

  idx1 <- which(all_x %in% lis1$x)
  idx2 <- which(all_x %in% lis2$x)
  idx_vec <- rep(NA, n)
  idx_vec[idx1] <- 1; idx_vec[idx2] <- 2

  y1 <- rep(0, n)
  y2 <- rep(0, n)

  for(i in 1:n){
    if(idx_vec[i] == 1){
      y1[i] = all_y[i]
      if(i > 1) y2[i] = y2[i-1]
    } else {
      y2[i] = all_y[i]
      if(i > 1) y1[i] = y1[i-1]
    }
  }

  max(abs(y1-y2))
}


