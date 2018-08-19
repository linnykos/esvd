est_func <- function(vec){
  p_est <- mean(vec)/var(vec)
  n_est <- mean(vec)*p_est/(1-p_est)

  list(size = n_est, prob = p_est)
}

# unknown r and p
x <- rnbinom(1000, 100, 0.9)
est_func(x)

# log series distribution
library(vcdExtra)

generate_data <- function(n, theta, lambda){
  #method 1
  alpha <- -1/log(1-theta)
  P <- theta/(1-theta); Q <- 1+P
  x1 <- rnbinom(n, size = alpha*lambda, prob = 1-P/Q)

  x2 <- sapply(1:n, function(x){
    M <- rpois(1, lambda)
    if(M == 0) return(0)
    Y <- vcdExtra::rlogseries(n = M, prob = theta)
    sum(Y)
  })

  list(x1 = x1, x2 = x2)
}


res <- generate_data(10000, 0.1, 20)

png(paste0("../figure/experiment/12_nbd_distribution.png"), height = 1200, width = 1200, res = 300, units = "px")
plot(sort(res$x1), sort(res$x2), asp = T, pch = 16,
     xlab = "CDF of NBD", ylab = "CDF of Poisson sum of LSD",
     col = rgb(0, 0, 0, 0.1))
lines(c(-1000,1000), c(-1000,1000), col = "red", lwd = 2, lty = 2)
graphics.off()

###################

#EM algorithm

e_step <- function(vec, theta, lambda){
  n <- length(vec)
  alpha <- -1/log(1-theta)
  z <- -((1-theta)/theta - alpha)

  m_vec <- sapply(1:n, function(i){
    alpha*lambda*sum(sapply(1:vec[i], function(k){
      1/(alpha*lambda + k - 1)
    }))
  })

  sum_z <- m_vec * z

  list(m_vec = m_vec, sum_z = sum_z)
}

m_step <- function(vec, m_vec, sum_z){
  n <- length(vec)
  theta <- sum(vec - m_vec)/sum(vec + sum_z - m_vec)
  lambda <- sum(m_vec)/n

  list(theta = theta, lambda = lambda)
}

em_alg <- function(vec, theta, lambda, tol = 1e-6, max_iter = 100){
  iter <- 1

  while(TRUE){
    old_theta <- theta; old_lambda <- lambda

    res <- e_step(vec, theta, lambda)
    res <- m_step(vec, res$m_vec, res$sum_z)

    theta <- min(res$theta, 1-tol); lambda <- res$lambda

    if((abs(theta - old_theta) < tol & abs(lambda - old_lambda) < tol)
       | iter > max_iter) break()
    # print(iter)

    iter <- iter + 1
  }

  alpha <- -1/log(1-theta)
  P <- theta/(1-theta); Q <- 1+P

  list(size = alpha*lambda, prob = 1-P/Q)
}

vec <- rnbinom(1000, 10, 0.5)
res <- em_alg(vec, theta = 0.5, lambda = 5)

####################

#simulation setup
n <- 100
r <- 5
p_vec <- exp(seq(log(0.001), log(0.4), length.out = 10))
trials <- 10

size_mat1 <- matrix(NA, ncol = length(p_vec), nrow = trials)
size_mat2 <- matrix(NA, ncol = length(p_vec), nrow = trials)
prob_mat1 <- matrix(NA, ncol = length(p_vec), nrow = trials)
prob_mat2 <- matrix(NA, ncol = length(p_vec), nrow = trials)

for(i in 1:length(p_vec)){
  cat('*')

  for(j in 1:trials){
    set.seed(j)
    vec <- rnbinom(n, size = r, prob = p_vec[i])

    res1 <- est_func(vec)
    res2 <- em_alg(vec, theta = 0.5, lambda = 5)

    size_mat1[j,i] <- res1$size
    size_mat2[j,i] <- res2$size
    prob_mat1[j,i] <- res1$prob
    prob_mat2[j,i] <- res2$prob
  }
}

png(paste0("../figure/experiment/12_nbd_comparison.png"), height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(NA, xlim = range(p_vec),
     ylim = range(c(size_mat1, size_mat2)),
     xlab = "Probability of failure",
     ylab = "Estimated size")
points(p_vec, colMeans(size_mat1),
       pch = 16, col = rgb(0.803, 0.156, 0.211))
arrows(p_vec, colMeans(size_mat1) - apply(size_mat1, 2, sd),
       p_vec, colMeans(size_mat1) + apply(size_mat1, 2, sd),
       length = 0.05, angle = 90, code = 3,
       col = rgb(0.803, 0.156, 0.211))

points(p_vec, colMeans(size_mat2),
       pch = 16, col = rgb(0.584, 0.858, 0.564))
arrows(p_vec, colMeans(size_mat2) - apply(size_mat2, 2, sd),
       p_vec, colMeans(size_mat2) + apply(size_mat2, 2, sd),
       length = 0.05, angle = 90, code = 3,
       col = rgb(0.584, 0.858, 0.564))

lines(c(-10,10), rep(5, 2), lwd = 2, lty = 2)

##

plot(NA, xlim = range(p_vec),
     ylim = range(c(prob_mat1, prob_mat2)),
     xlab = "Probability of failure",
     ylab = "Estimated probability")
points(p_vec, colMeans(prob_mat1),
       pch = 16, col = rgb(0.803, 0.156, 0.211))
arrows(p_vec, colMeans(prob_mat1) - apply(prob_mat1, 2, sd),
       p_vec, colMeans(prob_mat1) + apply(prob_mat1, 2, sd),
       length = 0.05, angle = 90, code = 3,
       col = rgb(0.803, 0.156, 0.211))

points(p_vec, colMeans(prob_mat2),
       pch = 16, col = rgb(0.584, 0.858, 0.564))
arrows(p_vec, colMeans(prob_mat2) - apply(prob_mat2, 2, sd),
       p_vec, colMeans(prob_mat2) + apply(prob_mat2, 2, sd),
       length = 0.05, angle = 90, code = 3,
       col = rgb(0.584, 0.858, 0.564))

lines(c(-10,10), c(-10, 10), lwd = 2, lty = 2)
graphics.off()

########
# other end of the spectrum

#simulation setup
n <- 100
r <- 5
p_vec <- 1-rev(exp(seq(log(0.001), log(0.4), length.out = 10)))
trials <- 10

size_mat1 <- matrix(NA, ncol = length(p_vec), nrow = trials)
prob_mat1 <- matrix(NA, ncol = length(p_vec), nrow = trials)

for(i in 1:length(p_vec)){
  cat('*')

  for(j in 1:trials){
    set.seed(j)
    vec <- rnbinom(n, size = r, prob = p_vec[i])

    res1 <- est_func(vec)
    size_mat1[j,i] <- res1$size
    prob_mat1[j,i] <- res1$prob
  }
}

p_vec <- p_vec[1:7]
size_mat1 <- size_mat1[,1:7]
prob_mat1 <- prob_mat1[,1:7]

png(paste0("../figure/experiment/12_nbd_comparison_2.png"), height = 1200, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2))
plot(NA, xlim = range(p_vec),
     ylim = range(size_mat1),
     xlab = "Probability of failure",
     ylab = "Estimated size")
points(p_vec, colMeans(size_mat1),
       pch = 16, col = rgb(0.803, 0.156, 0.211))
arrows(p_vec, colMeans(size_mat1) - apply(size_mat1, 2, sd),
       p_vec, colMeans(size_mat1) + apply(size_mat1, 2, sd),
       length = 0.05, angle = 90, code = 3,
       col = rgb(0.803, 0.156, 0.211))

lines(c(-10,10), rep(5, 2), lwd = 2, lty = 2)
lines(c(-10,10), rep(0, 2), lwd = 2)

##

plot(NA, xlim = range(p_vec),
     ylim = range(prob_mat1),
     xlab = "Probability of failure",
     ylab = "Estimated probability")
points(p_vec, colMeans(prob_mat1),
       pch = 16, col = rgb(0.803, 0.156, 0.211))
arrows(p_vec, colMeans(prob_mat1) - apply(prob_mat1, 2, sd),
       p_vec, colMeans(prob_mat1) + apply(prob_mat1, 2, sd),
       length = 0.05, angle = 90, code = 3,
       col = rgb(0.803, 0.156, 0.211))

lines(c(-10,10), c(-10, 10), lwd = 2, lty = 2)
lines(c(-10,10), rep(1, 2), lwd = 2)
graphics.off()

############################

# variance of the estimated mean and estimated variance
p_vec <- seq(0.01, .99, length.out = 20)
trials <- 50
n <- 100
r <- 5

mean_vec <- rep(NA, length(p_vec))
var_vec <- rep(NA, length(p_vec))

for(i in 1:length(p_vec)){
  tmp_mean <- rep(NA, trials)
  tmp_var <- rep(NA, trials)

  for(j in 1:trials){
    vec <- rnbinom(n, r, p_vec[i])
    tmp_mean[j] <- mean(vec)
    tmp_var[j] <- var(vec)
  }

  mean_vec[i] <- var(tmp_mean)
  var_vec[i] <- var(tmp_var)
}

