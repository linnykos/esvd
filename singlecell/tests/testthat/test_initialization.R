context("Test initialization")

## .svd_projection is correct

## .adaptive_gradient_step is correct

## .ensure_feasibility is correct

####################

## .projection_l1 is correct

test_that(".projection_l1 works", {
  set.seed(10)
  current_vec <- rnorm(10)
  other_mat <- matrix(rnorm(50), ncol = 10)
  res <- .projection_l1(current_vec, other_mat, other_bound = -100)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == length(current_vec))
})

test_that(".projection_l1 enforces negativity", {
  set.seed(10)
  current_vec <- rnorm(10)
  other_mat <- matrix(rnorm(50), ncol = 10)
  res <- .projection_l1(current_vec, other_mat, other_bound = -100)

  expect_true(all(other_mat %*% res <= 0))
})

test_that(".projection_l1 enforces positivity", {
  set.seed(10)
  current_vec <- rnorm(10)
  other_mat <- matrix(rnorm(50), ncol = 10)
  res <- .projection_l1(current_vec, other_mat, direction = ">=", other_bound = 100)

  expect_true(all(other_mat %*% res >= 0))
})

test_that(".projection_l1 maintains less than 0", {
  trials <- 50

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)

    dat <- abs(matrix(rexp(40, 1), nrow = 10, ncol = 4))

    res <- initialization(dat, max_val = -100)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    prod_mat <- u_mat %*% t(v_mat)

    all(prod_mat[which(!is.na(dat))] <= 0)
  })

  expect_true(all(bool_vec))
})

test_that(".projection_l1 maintains greater than 0", {
  trials <- 50

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)

    dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))

    res <- initialization(dat, family = "gaussian", max_val = 100)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    prod_mat <- u_mat %*% t(v_mat)

    all(prod_mat[which(!is.na(dat))] > 0)
  })

  expect_true(all(bool_vec))
})

test_that(".projection_l1 can keep the current vector", {
  set.seed(20)
  current_vec <- abs(rnorm(10))
  other_mat <- -abs(matrix(rnorm(50), ncol = 10))
  res <- .projection_l1(current_vec, other_mat, other_bound = -100)

  expect_true(sum(abs(res - current_vec)) < 1e-6)
})

test_that(".projection_l1 is actually a projection compared to the all 0 vector", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    current_vec <- abs(rnorm(10))
    other_mat <- -abs(matrix(rnorm(50), ncol = 10))
    new_vec <- .projection_l1(current_vec, other_mat, other_bound = -100)

    sum(abs(current_vec - new_vec)) <= sum(abs(current_vec))
  })

  expect_true(all(bool_vec))
})

test_that(".projection_l1 can take another bound", {
  set.seed(10)
  current_vec <- rnorm(10)
  other_mat <- matrix(rnorm(50), ncol = 10)
  res <- .projection_l1(current_vec, other_mat, direction = ">=", other_bound = 1)

  expect_true(all(other_mat %*% res <= 1+1e-6))
})

