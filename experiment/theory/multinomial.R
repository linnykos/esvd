trials <- 1000

vec <- sapply(1:trials, function(x){
  set.seed(x)
  tmp1 <- runif(1); tmp2 <- runif(1)
  a <- min(tmp1, tmp2)
  b <- max(tmp1, tmp2)

  a <= (1-b)/(1-a)
})

idx <- which(!vec)[1]

x <- 6
set.seed(x)
tmp1 <- runif(1); tmp2 <- runif(1)
a <- min(tmp1, tmp2)
b <- max(tmp1, tmp2)

a <= (1-b)/(1-a)
