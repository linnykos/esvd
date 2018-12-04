context("Test initialization")

## initialization is correct

test_that("initialization works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))

  res <- initialization(dat)

  expect_true(is.list(res))
  expect_true(ncol(res$u_mat) == ncol(res$v_mat))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
})

test_that("initialization works off of an imputed matrix", {
  set.seed(10)
  dat <- abs(rbind(MASS::mvrnorm(5, rep(0, 5), diag(5)),
                   MASS::mvrnorm(5, rep(10, 5), 2*diag(5))))
  dat2 <- dat
  for(i in 1:nrow(dat2)){
    dat2[i, sample(1:5, 1)] <- NA
  }
  drop_idx <- which(is.na(dat2))
  dat2 <- .scImpute(dat, drop_idx, Kcluster = 2, min_size = 3)

  res <- initialization(dat2)

  expect_true(is.list(res))
  expect_true(ncol(res$u_mat) == ncol(res$v_mat))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
})

test_that("initialization gives negative predictions", {
  set.seed(10)
  dat <- abs(rbind(MASS::mvrnorm(5, rep(0, 5), diag(5)),
                   MASS::mvrnorm(5, rep(10, 5), 2*diag(5))))
  dat2 <- dat
  for(i in 1:nrow(dat2)){
    dat2[i, sample(1:5, 1)] <- NA
  }
  drop_idx <- which(is.na(dat2))
  dat2 <- .scImpute(dat, drop_idx, Kcluster = 2, min_size = 3)

  res <- initialization(dat2)
  pred_mat <- res$u_mat %*% t(res$v_mat)

  expect_true(all(pred_mat[which(!is.na(dat))] <= 1e-6))
})

test_that("initialization can catch hard cases", {
  # here, I think the low-rank U_mat's column space does not include ANY
  ## point that are all positive...

  set.seed(10*2)
  dat <- abs(matrix(rpois(25, 2), nrow = 5, ncol = 5))
  class(dat) <- c("poisson", class(dat)[length(class(dat))])
  init <- initialization(dat, family = "poisson", max_val = 100)

  expect_true(length(init) == 2)
})

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

