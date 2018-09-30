rm(list=ls())
obj_func <- function(vec, param1, param2){
  param <- param1*param2
  sum(log(param) + 2 * vec^2 /param^2 - 4 * vec/param)
}

obj_func2 <- function(vec, param1, param2){
  param <- param1*param2
  val <- 0
  for(i in 1:length(vec)){
    val <- val + log(param) + 2*vec[i]^2/param^2 - 4*vec[i]/param
  }
  val
}

grad_func <- function(vec, param1, param2){
  param <- param1*param2
  sum(param2*(1/param - 4*vec^2/param^3 + 4*vec/param^2))
}

grad_func2 <- function(vec, param1, param2){
  param <- param1*param2
  val <- 0
  for(i in 1:length(vec)){
    val <- val + param2/param - 4*vec[i]^2*param2/param^3 + 4*vec[i]*param2/param^2
  }
  val
}


trials <- 100

bool_vec <- sapply(1:trials, function(x){
  set.seed(x)
  val <- abs(rnorm(1, mean = 10, sd = 10))
  vec <- matrix(rnorm(100, mean = val, sd = val/2))
  param1 <- abs(rnorm(1, mean = 5, sd = 5))
  param1b <- abs(rnorm(1, mean = 5, sd = 5))
  param2 <- abs(rnorm(1, mean = 5, sd = 5))

  grad <- grad_func(vec, param1, param2)

  res <- obj_func(vec, param1, param2)
  res2 <- obj_func(vec, param1b, param2)

  res2 >= res + grad * (param1b - param1) - 1e-6
})
