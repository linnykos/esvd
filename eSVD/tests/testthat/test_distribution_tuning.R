context("Test distribution tuning")

## .compute_variance is correct

test_that(".compute_variance works", {
  mean_mat <- matrix(1:6, 3, 2)
  scalar <- 2
  res <- .compute_variance(mean_mat, family = "neg_binom", scalar = scalar)

  expect_true(all(dim(res) == dim(mean_mat)))
})

test_that(".compute_variance formula is correct for negative binomial", {
  mean_mat <- matrix(1:6, 3, 2)
  scalar <- 2
  res <- .compute_variance(mean_mat, family = "neg_binom", scalar = scalar)

  expect_true(all(res == mean_mat + mean_mat^2/scalar))
})

test_that(".compute_variance formula is correct for curved gaussian", {
  mean_mat <- matrix(1:6, 3, 2)
  scalar <- 2
  res <- .compute_variance(mean_mat, family = "curved_gaussian", scalar = scalar)

  expect_true(all(res == mean_mat^2/scalar^2))
})

###########

## .recompute_mean is correct

test_that(".recompute_mean works", {
  nat_mat <- -matrix(seq(0.1, 1, length.out = 10), 5, 2)
  mean_mat <- compute_mean(nat_mat, family = "exponential")
  scalar <- 1
  res <- .recompute_mean(nat_mat, mean_mat, scalar, recompute_mean = F)

  expect_true(all(res == mean_mat))
})

test_that(".recompute_mean will crash if any is NA", {
  nat_mat <- -matrix(seq(0.1, 1, length.out = 10), 5, 2)
  mean_mat <- compute_mean(nat_mat, family = "exponential")
  scalar <- 1
  mean_mat[1,1] <- NA
  expect_error(.recompute_mean(nat_mat, mean_mat, scalar, recompute_mean = F))
})

test_that(".recompute_mean can recompute mean", {
  nat_mat <- -matrix(seq(0.1, 1, length.out = 10), 5, 2)
  scalar <- 1
  mean_mat <- compute_mean(nat_mat, family = "neg_binom", scalar = scalar)
  res <- .recompute_mean(nat_mat, mean_mat, scalar, recompute_mean = T)

  expect_true(all(res == mean_mat))
})

test_that(".recompute_mean respects different scalars", {
  nat_mat <- -matrix(seq(0.1, 1, length.out = 10), 5, 2)
  mean_mat <- compute_mean(nat_mat, family = "neg_binom", scalar = 10)
  res <- .recompute_mean(nat_mat, mean_mat, scalar = 20, recompute_mean = T)

  expect_true(sum(abs(mean_mat - res)) > 0)
})

#########################

## .compute_tuning_objective is correct

test_that(".compute_tuning_objective works", {
  set.seed(10)
  scalar <- 50
  prob <- 0.25
  dat <- matrix(stats::rnbinom(30, size = scalar, prob = 1-prob), 5, 6)
  family <- "neg_binom"
  nat_mat <- matrix(log(prob), 5, 6)
  mean_mat <- matrix(scalar*prob/(1-prob), 5, 6)

  res <- .compute_tuning_objective(dat, family, nat_mat, mean_mat, scalar, recompute_mean = F)

  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
  expect_true(res > 0)
})

test_that(".compute_tuning_objective can work with recompute_mean", {
  set.seed(10)
  scalar <- 50
  prob <- 0.25
  dat <- matrix(stats::rnbinom(30, size = scalar, prob = 1-prob), 5, 6)
  family <- "neg_binom"
  nat_mat <- matrix(log(prob), 5, 6)
  mean_mat <- matrix(scalar*prob/(1-prob), 5, 6)

  res <- .compute_tuning_objective(dat, family, nat_mat, mean_mat, scalar, recompute_mean = T)

  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
  expect_true(res > 0)
})
