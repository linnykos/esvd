context("Matrix factorization")

## .subgradient_row is correct

test_that(".subgradient_row works", {
  set.seed(10)
  dat <- matrix(rnorm(200), 20, 10)
  initial_vec <- rnorm(5)
  latent_mat <- matrix(rnorm(50), 10, 5)

  fixed_idx <- 7
  index_in <- c(1,3,4,6)
  index_out <- c(2,8,10)

  res <- .subgradient_vec(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out)

  expect_true(is.numeric(res))
  expect_true(length(res) == 5)
})

test_that(".subgradient_row gives a correct subgradient for one instance", {
  set.seed(10)
  dat <- matrix(rnorm(200), 20, 10)
  initial_vec <- rnorm(5)
  latent_mat <- matrix(rnorm(50), 10, 5)

  fixed_idx <- 7
  index_in <- c(1,3,4,6)
  index_out <- c(2,8,10)

  res <- .subgradient_vec(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out)

  bool_vec <- sapply(1:100, function(x){
    set.seed(11*x)
    y <- rnorm(5)

    f1 <- .evaluate_objective_single(dat, y, latent_mat, fixed_idx, index_in, index_out)
    f2 <- .evaluate_objective_single(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out)

    f1 >= f2 + res %*% (y-initial_vec)
  })

  expect_true(all(bool_vec))
})

###########

## .evaluate_objective_single is correct

test_that(".evaluate_objective_single works", {
  set.seed(10)
  dat <- matrix(rnorm(200), 20, 10)
  initial_vec <- rnorm(5)
  latent_mat <- matrix(rnorm(50), 10, 5)

  fixed_idx <- 7
  index_in <- c(1,3,4,6)
  index_out <- c(2,8,10)

  res <- .evaluate_objective_single(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out)

  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
  expect_true(res >= 0)
})

test_that(".evaluate_objective_single evaluates correctly", {
  set.seed(20)
  dat <- matrix(rnorm(200), 20, 10)
  initial_vec <- rnorm(5)
  latent_mat <- matrix(rnorm(50), 10, 5)

  fixed_idx <- 7
  index_in <- c(1,3,4,6)
  index_out <- c(2,8,10)

  # manual calculation
  first_term <- 0
  for(j in index_in){
    first_term <- first_term + (dat[fixed_idx, j] - initial_vec %*% latent_mat[j,])^2
  }
  second_term <- 0
  for(j in index_out){
    second_term <- second_term + max(0, initial_vec %*% latent_mat[j,])^2
  }
  res2 <- first_term/(2*length(index_in)) + second_term/(2*length(index_out))

  res <- .evaluate_objective_single(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out)

  expect_true(abs(res - res2) <= 1e-6)
})
