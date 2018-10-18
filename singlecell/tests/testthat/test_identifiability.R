context("Test identifiability")

test_that(".identification works", {
  set.seed(10)
  X <- MASS::mvrnorm(n = 100, mu = rep(0,5), Sigma = diag(5))
  Y <- MASS::mvrnorm(n = 100, mu = rep(1,5), Sigma = 2*diag(5))
  res <- .identification(X, Y)

  expect_true(is.list(res))
  expect_true(all(sapply(res, is.matrix)))
  expect_true(all(dim(res$X) == dim(X)))
  expect_true(all(dim(res$Y) == dim(Y)))
})

test_that(".identification is correct", {
  set.seed(20)
  X <- MASS::mvrnorm(n = 100, mu = rep(0,5), Sigma = diag(5))
  Y <- MASS::mvrnorm(n = 100, mu = rep(1,5), Sigma = 2*diag(5))
  res <- .identification(X, Y)

  cov_x <- t(res$X) %*% res$X
  cov_y <- t(res$Y) %*% res$Y

  expect_true(sum(abs(cov_x - cov_y)) <= 1e-6)
})
