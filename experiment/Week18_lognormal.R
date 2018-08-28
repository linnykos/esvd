trials <- 1000
mu <- -5
sig <- 1

mat <- sapply(ceiling(exp(seq(log(10), log(5000), length.out = 10))), function(n){
  res <- sapply(1:trials, function(x) {
    set.seed(10*x)
    y <- rnorm(n, mean = mu, sd = sig)
    z <- exp(y)

    target <- exp(mu + sig^2/2)
    res1 <- mean(z)
    res2 <- exp(mean(y) + var(y)/2)

    c(abs(res1 - target), abs(res2 - target))
  })

  apply(res, 1, mean)
})

mat


