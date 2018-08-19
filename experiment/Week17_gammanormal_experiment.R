source("../experiment/em_gamma_normal.R")

generate_data <- function(n, shape = 1, rate, mean, sd, prop = 0.5,
                          positive = T, min_val = log10(1.01)){
  assign_vec <- sample(2, n, replace = T, prob = c(prop, 1-prop))
  x <- numeric(n)

  idx1 <- which(assign_vec == 1)
  x[idx1] <- stats::rgamma(length(idx1), shape = shape, rate = rate)

  idx2 <- which(assign_vec == 2)
  x[idx2] <- stats::rnorm(length(idx2), mean = mean, sd = sd)

  if(positive){
    x[x < min_val] <- min_val
  }

  x
}

x <- generate_data(1000, 1, .5, 4, 1, 0.5)
res <- get_mix(x)
