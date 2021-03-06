context("Test identifiability")

test_that(".identification works", {
  res <- .identification(diag(5), 2*diag(5))

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == 5))
})

test_that(".identification works on an instance where you need to worry about underflow", {
  load("../assets/identification1.RData")

  res <- .identification(cov_x, cov_y)
  expect_true(sum(abs(Im(res))) == 0)
})

test_that(".identification is correct for a deterministic setting", {
  cov_x <- matrix(c(2,1,1,2),2,2)
  cov_y <- matrix(c(5,-1,-1,5),2,2)
  res <- .identification(cov_x, cov_y)

  mat1 <- res %*% cov_x %*% t(res)
  mat2 <- solve(t(res)) %*% cov_y %*% solve(res)

  res <- .identification(diag(5), 2*diag(5))

  expect_true(sum(abs(mat1 - mat2)) <= 1e-6)
})

test_that(".identification issues warning for rank defficient setting", {
  set.seed(10)
  u_mat <- matrix(rnorm(12), ncol = 3, nrow = 4)
  v_mat <- matrix(rnorm(15), ncol = 3, nrow = 5)

  cov_x <- stats::cov(u_mat)
  cov_y <- stats::cov(v_mat)

  cov_y[,3] <- 0; cov_y[3,] <- 0
  expect_warning(.identification(cov_x, cov_y))
})

test_that(".identification is correct", {
  set.seed(20)
  cov_x <- cov(MASS::mvrnorm(n = 10, rep(0, 5), diag(5)))
  cov_y <- cov(MASS::mvrnorm(n = 10, rep(0, 5), toeplitz(5:1)))

  res <- .identification(cov_x, cov_y)

  cov_x2 <- res %*% cov_x %*% t(res)
  cov_y2 <- solve(t(res)) %*% cov_y %*% solve(res)

  expect_true(sum(abs(cov_x2 - cov_y2)) <= 1e-6)
  expect_true(abs(sum(cov_x2) - sum(diag(cov_x2))) <= 1e-6)
  expect_true(abs(sum(cov_y2) - sum(diag(cov_y2))) <= 1e-6)
})

test_that(".identification for 1-dim covariances (just variances)", {
  cov_x <- matrix(10, 1, 1)
  cov_y <- matrix(5, 1, 1)

  res <- .identification(cov_x, cov_y)

  cov_x2 <- res %*% cov_x %*% t(res)
  cov_y2 <- solve(t(res)) %*% cov_y %*% solve(res)

  expect_true(sum(abs(cov_x2 - cov_y2)) <= 1e-6)
})

####################

## .reparameterize is correct

test_that(".reparameterize works", {
  set.seed(10)
  u_mat <- MASS::mvrnorm(60, rep(0, 5), diag(5))
  v_mat <- MASS::mvrnorm(50, rep(1, 5), 2*diag(5))

  res <- .reparameterize(u_mat, v_mat)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(dim(res$u_mat) == dim(u_mat)))
  expect_true(all(dim(res$v_mat) == dim(v_mat)))
})

test_that(".reparameterize works for rank 1 matrices", {
  set.seed(10)
  u_mat <- matrix(rnorm(50), ncol = 1)
  v_mat <-matrix(rnorm(50), ncol = 1)

  res <- .reparameterize(u_mat, v_mat)

  expect_true(length(res) == 2)
})

test_that(".reparameterize preserves the inner products", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    u_mat <- MASS::mvrnorm(10, rep(0, 5), diag(5))
    v_mat <- MASS::mvrnorm(10, rep(1, 5), 2*diag(5))

    res <- .reparameterize(u_mat, v_mat)

    pred_mat <- u_mat %*% t(v_mat)
    pred_mat2 <- res$u_mat %*% t(res$v_mat)

    sum(abs(pred_mat - pred_mat2)) <= 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".reparameterize yields the same second moment matrix", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(11*x)
    u_mat <- MASS::mvrnorm(5, rep(0, 5), diag(5))
    v_mat <- MASS::mvrnorm(5, rep(1, 5), 2*diag(5))

    res <- .reparameterize(u_mat, v_mat)

    cov_u <- t(res$u_mat) %*% res$u_mat
    cov_v <- t(res$v_mat) %*% res$v_mat

    sum(abs(cov_u - cov_v)) < 1e-6
  })

  expect_true(all(bool_vec))
})


test_that(".reparameterize yields the same covariance matrix", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(12*x)
    n <- 10
    p <- 20
    u_mat <- MASS::mvrnorm(n, rep(0, 5), diag(5))
    v_mat <- MASS::mvrnorm(p, rep(1, 5), 2*diag(5))

    res <- .reparameterize(u_mat, v_mat, equal_covariance = T)

    cov_u <- t(res$u_mat) %*% res$u_mat/n
    cov_v <- t(res$v_mat) %*% res$v_mat/p

    sum(abs(cov_u - cov_v)) < 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".reparameterize yields diagonal covariances", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(13*x)
    u_mat <- MASS::mvrnorm(10, rep(0, 5), diag(5))
    v_mat <- MASS::mvrnorm(20, rep(1, 5), 2*diag(5))

    res <- .reparameterize(u_mat, v_mat)

    cov_u <- t(res$u_mat) %*% res$u_mat
    cov_v <- t(res$v_mat) %*% res$v_mat

    max(abs(sum(cov_u) - sum(diag(cov_u))), abs(sum(cov_v) - sum(diag(cov_v)))) < 1e-6
  })

  expect_true(all(bool_vec))
})

