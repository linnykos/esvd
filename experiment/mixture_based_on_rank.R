rm(list = ls())
.random_walk <- function(k){
  x <- sample(c(1, -1), k, replace = T)
  #x <- sign(cumsum(x))
}

.l2norm <- function(x){
  sqrt(sum(x^2))
}

.generate_latent <- function(n = 100, k = 2){
  dat <- t(sapply(1:n, function(x){
    .random_walk(k)
  }))

  if(k == 1){
    dat <- matrix(dat, ncol = 1)
  }

  res <- t(apply(dat, 1, function(x){
    if(k == 1){
      MASS::mvrnorm(n = 1, mu = as.matrix(x), Sigma = as.matrix(.25))
    } else {
      MASS::mvrnorm(n = 1, mu = x, Sigma = .1*diag(k))
    }
  }))

  if(k == 1){
    res <- matrix(res, ncol = 1)
  }

  for(i in 1:nrow(res)){
    res[i,] <- res[i,]/.l2norm(res[i,])
  }

  res
}

.generate_samples <- function(mat){
  gram <- mat %*% t(mat)

  vec <- gram[lower.tri(gram, diag = F)]
  sapply(vec, function(x){
    stats::rnorm(1, mean = x, sd = .05)
  })
}

set.seed(10)
latent <- .generate_latent(k = 4)
dat <- .generate_samples(latent)
hist(dat, breaks = 50, col = "gray")
