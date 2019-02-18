context("Test identifiability")

test_that(".identification works", {
  res <- .identification(diag(5), 2*diag(5))

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == 5))
})

test_that(".identification is correct", {
  set.seed(20)
  cov_x <- cov(MASS::mvrnorm(n = 10, rep(0, 5), diag(5)))
  cov_y <- cov(MASS::mvrnorm(n = 10, rep(0, 5), toeplitz(5:1)))

  res <- .identification(cov_x, cov_y)

  cov_x2 <- res %*% cov_x %*% t(res)
  cov_y2 <- solve(t(res)) %*% cov_y %*% solve(res)

  expect_true(sum(abs(cov_x2 - cov_y2)) <= 1e-6)
})
