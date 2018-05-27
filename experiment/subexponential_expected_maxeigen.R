rm(list=ls())
set.seed(10)

#fix a vector
n <- 100
trials <- 1000
vec <- rep(0, trials)
vec2 <- rep(0, trials)

v <- rnorm(n)
v <- v/sqrt(sum(v^2))

for(k in 1:trials){
  if(k %% floor(trials/10) == 0) cat('*')

  set.seed(k)
  X <- matrix(0, n, n)
  for(i in 1:n){
    for(j in i:n){
      X[i,j] <- rexp(1, rate = 1/2) - 2
      X[j,i] <- X[i,j]
    }
  }

  vec[k] <- t(v)%*%X%*%v
  vec2[k] <- max(eigen(X)$values)
}

hist(vec, breaks = 50, col = "gray")

#############

# test #2
n_vec <- exp(seq(log(10), log(1000), length.out = 10))
trials <- 1000
vec_prob <- rep(0, length(n_vec))
vec_exp <- rep(0, length(n_vec))

for(n in 1:length(n_vec)){
  print('*')
  tmp <- rep(0, trials)
  for(k in 1:trials){
    set.seed(k)
    X <- matrix(0, n_vec[n], n_vec[n])
    for(i in 1:n){
      for(j in i:n){
        X[i,j] <- rexp(1, rate = 1/2) - 2
        X[j,i] <- X[i,j]
      }
    }

    tmp[k] <- max(eigen(X)$values)
  }

  vec_prob[n] <- quantile(tmp, probs = 1 - 1/n_vec[n])
  vec_exp[n] <- mean(tmp)
}

plot(n_vec, vec_prob)
zz <- sqrt(n_vec + log(n_vec))
const <- lm(vec_prob ~ zz)
const <- stats::coef(const)
lines(n_vec, const[1]+const[2]* zz, col = "red", lwd = 2)
