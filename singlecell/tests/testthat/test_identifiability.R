context("Test identifiability")

test_that(".identification works", {
  res <- .identification(diag(5), 2*diag(5))

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == 5))
})

test_that(".identification is correct", {
  set.seed(20)
  cov1 <- cov(MASS::mvrnorm(n = 10, rep(0, 5), diag(5)))
  cov2 <- cov(MASS::mvrnorm(n = 10, rep(0, 5), 2*diag(5)))

  res <- .identification(cov1, cov2)

  cov_x <- res %*% cov1 %*% t(res)
  cov_y <- solve(t(res)) %*% cov2 %*% solve(res)

  expect_true(sum(abs(cov_x - cov_y)) <= 1e-6)
})

test_that(".identification preserves the inner product", {
  set.seed(20)
  X <- MASS::mvrnorm(n = 100, mu = rep(0,5), Sigma = diag(5))
  Y <- MASS::mvrnorm(n = 100, mu = rep(1,5), Sigma = 2*diag(5))
  res <- .identification(X, Y)

  pred_mat <- X %*% t(Y)
  pred_mat2 <- res$X %*% t(res$Y)

  expect_true(sum(abs(pred_mat - pred_mat2)) <= 1e-6)
})
