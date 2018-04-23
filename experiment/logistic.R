trials <- 1000

dat <- sapply(1:trials, function(x){
  set.seed(x)
  vec1 <- rnorm(5)
  vec2 <- rnorm(5)
  val <- vec1 %*% vec2

  1/(1+exp(-val))
})

plot(sort(dat))
