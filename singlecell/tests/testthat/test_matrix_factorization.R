context("Test matrix factorization")

# .gradient_vec is correct

test_that(".gradient_vec works", {
  set.seed(10)
  dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)

  i <- 5
  res <- .gradient_vec(dat[i,], u_mat[i,], v_mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == 2)
})

test_that(".gradient_vec works for the other direction", {
  set.seed(8)
  dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)

  j <- 2
  res <- .gradient_vec(dat[,j], v_mat[j,], u_mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == 2)
})

##################

## .evaluate_objective is correct

test_that(".evaluate_objective works", {
  set.seed(20)
  dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)

  res <- .evaluate_objective(dat, u_mat, v_mat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(!is.nan(res))
})

test_that(".evaluate_objective yields a smaller value under truth", {
  set.seed(20)
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)
  dat <- u_mat %*% t(v_mat)
  dat[sample(1:prod(dim(dat)), 10)] <- NA

  u_mat2 <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat2 <- matrix(rnorm(8), nrow = 4, ncol = 2)

  res <- .evaluate_objective(dat, u_mat, v_mat)
  res2 <- .evaluate_objective(dat, u_mat2, v_mat2)

  expect_true(res < res2)
})

