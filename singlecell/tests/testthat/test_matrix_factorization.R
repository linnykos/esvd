context("Test matrix factorization")

# .gradient_vec is correct

test_that(".gradient_vec works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
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
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)

  j <- 2
  res <- .gradient_vec(dat[,j], v_mat[j,], u_mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == 2)
})

test_that(".gradient_vec direction lowers objective value", {
  set.seed(15)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)

  i <- 5
  u_vec <- u_mat[i,]
  grad <- .gradient_vec(dat[i,], u_vec, v_mat)
  u_vec2 <- u_vec - 0.01*grad

  res <- .evaluate_objective_single(dat[i,], u_vec, v_mat)
  res2 <- .evaluate_objective_single(dat[i,], u_vec2, v_mat)

  expect_true(res2 < res)
})

test_that(".gradient_vec satisfies the gradient definition", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
    u_vec <- abs(rnorm(2)); u_vec2 <- abs(rnorm(2))
    v_mat <- -abs(matrix(rnorm(8), nrow = 4, ncol = 2))

    i <- sample(1:10, 1)
    grad <- .gradient_vec(dat[i,], u_vec, v_mat)

    res <- .evaluate_objective_single(dat[i,], u_vec, v_mat)
    res2 <- .evaluate_objective_single(dat[i,], u_vec2, v_mat)

    res2 >= res + as.numeric(grad %*% (u_vec2 - u_vec))
  })

  expect_true(all(bool_vec))
})


##################

## .evaluate_objective is correct

test_that(".evaluate_objective works", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- -matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)

  res <- .evaluate_objective(dat, u_mat, v_mat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(!is.nan(res))
})

test_that(".evaluate_objective yields a smaller value under truth", {
  set.seed(20)
  dat <- matrix(1, nrow = 10, ncol = 4)
  svd_res <- svd(dat)
  u_mat <- -svd_res$u %*% diag(sqrt(svd_res$d))
  v_mat <- svd_res$v %*% diag(sqrt(svd_res$d))
  dat[sample(1:prod(dim(dat)), 10)] <- NA

  u_mat2 <- -matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat2 <- matrix(rnorm(8), nrow = 4, ncol = 2)

  res <- .evaluate_objective(dat, u_mat, v_mat)
  res2 <- .evaluate_objective(dat, u_mat2, v_mat2)

  expect_true(res < res2)
})

################

## .evaluate_objective_single is correct

test_that(".evaluate_objective_single works", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- -matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)

  i <- 5
  res <- .evaluate_objective_single(dat[i,], u_mat[i,], v_mat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(!is.nan(res))
})

test_that(".evaluate_objective_single yields a smaller value under truth", {
  set.seed(20)
  dat <- matrix(1, nrow = 10, ncol = 4)
  svd_res <- svd(dat)
  u_mat <- -svd_res$u %*% diag(sqrt(svd_res$d))
  v_mat <- svd_res$v %*% diag(sqrt(svd_res$d))
  dat[sample(1:prod(dim(dat)), 10)] <- NA

  u_mat2 <- matrix(rnorm(20), nrow = 10, ncol = 4)

  i <- 5
  res <- .evaluate_objective_single(dat[i,], u_mat[i,], v_mat)
  res2 <- .evaluate_objective_single(dat[i,], u_mat2[i,], v_mat)

  expect_true(res < res2)
})

