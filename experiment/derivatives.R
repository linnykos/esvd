rm(list=ls())
x <- seq(-0.1, -20, length.out = 100)
alpha <- 0.1

y <- -log(-x) + alpha*x + alpha^2*x^2

plot(x, y, asp = T)

# plot the derivative
z <- -1/x + alpha + 2*alpha^2*x
plot(x, z, asp = T)

# plot the second derivative
s <- 1/x^2 + 2*alpha^2
plot(x, s, asp = T, ylim = c(-10,10))

######################

set.seed(10)
trials <- 10000
mu <- 5
alpha <- 4
x <- rnorm(trials, mean = mu, sd = mu/alpha)
var(x)
(mu/alpha)^2
x2 <- x^2
var(x2)
2*(1+2*alpha^2)*((mu/alpha)^4)
# https://stats.stackexchange.com/questions/119818/non-zero-mean-and-finite-variance-gaussian-squared-r-v-has-non-central-chi-squar

