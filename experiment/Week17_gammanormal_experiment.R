source("../experiment/em_gamma_normal.R")
source("../experiment/em_gamma_truncatednormal.R")
source("../experiment/truncated_normal_estimation.R")

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

set.seed(10)
x <- generate_data(1000, 1, .5, 0.5, 4, 1)
res_original <- get_mix(x)
res_trunc <- .get_mix(x)


###################

vec1 = rnorm(100); vec1 = vec1[vec1 > 0]
vec2 = rnorm(100); vec2 = vec2[vec2 > 0]
vec3 = rnorm(100, mean = 4); vec3 = vec3[vec3 > 0]

res1 = .estimate_cdf(vec1)
res2 = .estimate_cdf(vec2)
res3 = .estimate_cdf(vec3)

.ks_distance(res1, res2)
.ks_distance(res1, res3)
