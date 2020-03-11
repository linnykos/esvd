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

test_that(".recompute_mean computes the mean correctly", {
  nat_mat <- -matrix(seq(0.1, 1, length.out = 10), 5, 2)
  scalar <- 10
  mean_mat <- compute_mean(nat_mat, family = "neg_binom", scalar = 2)
  res <- .recompute_mean(nat_mat, mean_mat, scalar, recompute_mean = T)

  expect_true(sum(abs(res - scalar*exp(nat_mat)/(1-exp(nat_mat)))) <= 1e-6)
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

test_that(".compute_tuning_objective is a meaningful ratio when scalar is correct", {
  set.seed(10)
  scalar <- 50
  prob <- 0.25
  dat <- matrix(stats::rnbinom(200, size = scalar, prob = 1-prob), 10, 20)
  nat_mat <- matrix(log(prob), nrow = nrow(dat), ncol = ncol(dat))
  mean_mat <- compute_mean(nat_mat, family = "neg_binom", scalar = scalar)

  res <- .compute_tuning_objective(dat, family = "neg_binom", nat_mat = nat_mat,
                                   mean_mat = mean_mat, scalar = scalar, recompute_mean = F)

  res2 <- .compute_tuning_objective(dat, family = "neg_binom", nat_mat = nat_mat,
                                   mean_mat = mean_mat, scalar = 1, recompute_mean = F)

  res3 <- .compute_tuning_objective(dat, family = "neg_binom", nat_mat = nat_mat,
                                    mean_mat = mean_mat, scalar = 2000, recompute_mean = F)

  expect_true(abs(res - prod(dim(dat))) < abs(res2 - prod(dim(dat))))
  expect_true(abs(res - prod(dim(dat))) < abs(res3 - prod(dim(dat))))
})

##################

## .tuning_fit is correct

test_that(".tuning_fit works", {
  set.seed(10)
  scalar <- 50
  prob <- 0.25
  dat <- matrix(stats::rnbinom(30, size = scalar, prob = 1-prob), 5, 6)

  res <- .tuning_fit(dat, family = "neg_binom", scalar = scalar,
                     max_val = 100, max_iter = 10, k = 2)

  expect_true(is.list(res))
  expect_true(all(c("u_mat", "v_mat") %in% names(res)))
})

test_that(".tuning_fit can work under randomly generated rank 1 cases", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    u_mat <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
    v_mat <- -abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    pred_mat <- u_mat %*% t(v_mat)
    dat <- pred_mat

    for(i in 1:10){
      for(j in 1:4){
        dat[i,j] <- stats::rnbinom(1, size = 10, prob = 1-exp(pred_mat[i,j]))
      }
    }

    res <- .tuning_fit(dat, family = "neg_binom", scalar = 10, max_val = 100,
                       max_iter = 10, k = 2)

    is.list(res) & c("u_mat", "v_mat") %in% names(res)
  })

  expect_true(all(bool_vec))
})

##############################

## .tuning_param_search is correct

test_that(".tuning_param_search works", {
  set.seed(10)
  scalar <- 50
  prob <- 0.25
  dat <- matrix(stats::rnbinom(200, size = scalar, prob = 1-prob), 10, 20)
  u_mat <- matrix(1, nrow = nrow(dat), ncol = 1)
  v_mat <- matrix(log(prob), nrow = ncol(dat), ncol = 1)

  res <- .tuning_param_search(dat, u_mat, v_mat, family = "neg_binom")

  expect_true(res > 1)
})

test_that(".tuning_param_search works on a more realistic setting for poisson", {
  set.seed(10)
  u_mat <- abs(matrix(rnorm(60), nrow = 30, ncol = 2))
  v_mat <- -abs(matrix(rnorm(60), nrow = 30, ncol = 2))
  pred_mat <- u_mat %*% t(v_mat)
  dat <- pred_mat

  for(i in 1:30){
    for(j in 1:30){
      dat[i,j] <- stats::rnbinom(1, size = 50, prob = 1-exp(pred_mat[i,j]))
    }
  }

  init <- .tuning_fit(dat, family = "poisson", max_val = 100, max_iter = 10, k = 1)
  res <- .tuning_param_search(dat, init$u_mat, init$v_mat, family = "poisson")

  expect_true(res > 1)
})

test_that(".tuning_param_search works on a more realistic setting for negative binomial", {
  set.seed(10)
  u_mat <- abs(matrix(rnorm(60), nrow = 30, ncol = 2))
  v_mat <- -abs(matrix(rnorm(60), nrow = 30, ncol = 2))
  pred_mat <- u_mat %*% t(v_mat)
  dat <- pred_mat

  for(i in 1:30){
    for(j in 1:30){
      dat[i,j] <- stats::rnbinom(1, size = 50, prob = 1-exp(pred_mat[i,j]))
    }
  }

  init <- .tuning_fit(dat, family = "neg_binom", max_val = 100, max_iter = 10, k = 1, scalar = 50)
  res <- .tuning_param_search(dat, init$u_mat, init$v_mat, family = "neg_binom")

  expect_true(res > 1)
})

test_that(".tuning_param_search works on a more realistic setting for curved gaussian", {
  set.seed(10)
  u_mat <- abs(matrix(rnorm(60), nrow = 30, ncol = 2))
  v_mat <- abs(matrix(rnorm(60), nrow = 30, ncol = 2))
  pred_mat <- u_mat %*% t(v_mat)
  dat <- pred_mat

  for(i in 1:30){
    for(j in 1:30){
      dat[i,j] <- max(stats::rnorm(1, mean = 1/pred_mat[i,j], sd = 1/(2*pred_mat[i,j])), 1e-3)
    }
  }

  init <- .tuning_fit(dat, family = "curved_gaussian", max_val = 100, max_iter = 10, k = 1, scalar = 2)
  res <- .tuning_param_search(dat, init$u_mat, init$v_mat, family = "curved_gaussian")

  expect_true(res > 1)
})

#########################

## tuning_scalar works

test_that("tuning_scalar works", {
  set.seed(10)
  scalar <- 50
  prob <- 0.25
  dat <- matrix(stats::rnbinom(200, size = scalar, prob = 1-prob), 10, 20)
  res <- tuning_scalar(dat, family = "neg_binom", max_val = 100, k = 1)

  expect_true(is.numeric(res))
  expect_true(all(res >= 1))
})

test_that("tuning_scalar works for rank 1 negative binomial", {
  set.seed(10)
  u_mat <- abs(matrix(rnorm(60), nrow = 30, ncol = 2))
  v_mat <- -abs(matrix(rnorm(60), nrow = 30, ncol = 2))
  pred_mat <- u_mat %*% t(v_mat)
  dat <- pred_mat

  for(i in 1:30){
    for(j in 1:30){
      dat[i,j] <- stats::rnbinom(1, size = 50, prob = 1-exp(pred_mat[i,j]))
    }
  }

  res <- tuning_scalar(dat, family = "neg_binom", max_val = 100, k = 1)

  expect_true(is.numeric(res))
  expect_true(all(res >= 1))
})

test_that("tuning_scalar works for rank 1 curved gaussian", {
  set.seed(10)
  u_mat <- abs(matrix(rnorm(60), nrow = 30, ncol = 2))
  v_mat <- abs(matrix(rnorm(60), nrow = 30, ncol = 2))
  pred_mat <- u_mat %*% t(v_mat)
  dat <- pred_mat

  for(i in 1:30){
    for(j in 1:30){
      dat[i,j] <- max(stats::rnorm(1, mean = 1/pred_mat[i,j], sd = 1/(2*pred_mat[i,j])), 1e-3)
    }
  }

  res <- tuning_scalar(dat, family = "curved_gaussian", max_val = 100, k = 1)

  expect_true(is.numeric(res))
  expect_true(all(res >= 1))
})
