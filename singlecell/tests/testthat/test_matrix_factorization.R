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

test_that(".subgradient_row for more settings", {
  bool_vec <- sapply(1:200, function(x){
    dat <- matrix(rnorm(200), 20, 10)
    initial_vec <- rnorm(5)
    latent_mat <- matrix(rnorm(50), 10, 5)

    fixed_idx <- 7
    index_in <- c(1,3,4,6)
    index_out <- c(2,8,10)

    res <- .subgradient_vec(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out)

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

#################

test_that(".evaluate_objective_single works", {
  set.seed(10)
  dat <- matrix(rnorm(200), 20, 10)
  initial_vec <- rnorm(5)
  latent_mat <- matrix(rnorm(50), 10, 5)

  fixed_idx <- 7
  index_in <- c(1,3,4,6)
  index_out <- c(2,8,10)

  res <- .estimate_row(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out)

  expect_true(length(res) == 5)
  expect_true(is.numeric(res))
})

test_that(".evaluate_objective_single works with verbose", {
  set.seed(10)
  dat <- matrix(rnorm(200), 20, 10)
  initial_vec <- rnorm(5)
  latent_mat <- matrix(rnorm(50), 10, 5)

  fixed_idx <- 7
  index_in <- c(1,3,4,6)
  index_out <- c(2,8,10)

  res <- .estimate_row(dat, initial_vec, latent_mat, fixed_idx, index_in, index_out,
                       max_iter = 100, verbose = T)

  expect_true(is.list(res))
  expect_true(length(res$obj_vec) == 101)
  expect_true(res$obj_vec[1] > res$obj_vec[101])
})

################

## .estimate_matrix is correct

test_that(".estimate_matrix works", {
  set.seed(10)
  dat <- matrix(rnorm(200), 20, 10)
  initial_mat <- matrix(rnorm(50), 20, 5)
  latent_mat <- matrix(rnorm(50), 10, 5)

  pattern <- matrix(0, 20, 10)
  pattern[sample(1:200, 100)] <- 1
  index_in_mat <- which(pattern == 1, arr.ind = T)
  index_out_mat <- which(pattern == 0, arr.ind = T)

  res <- .estimate_matrix(dat, initial_mat, latent_mat, index_in_mat, index_out_mat,
                          max_iter = 50)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == dim(initial_mat)))
})

#######################

## .convert_index_to_position is correct

test_that(".convert_index_to_position works", {
  vec <- 1:50
  num_row <- 10
  num_col <- 10

  res <- .convert_index_to_position(vec, num_row, num_col)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == c(length(vec), 2)))
  expect_true(all(res[,1] <= num_row))
  expect_true(all(res[,2] <= num_col))
})

test_that(".convert_index_to_position gives the same indices as which", {
  set.seed(10)
  mat <- matrix(0, 20, 10)
  mat[sample(1:200, 100)] <- 1

  vec1 <- which(mat == 1)
  res <- .convert_index_to_position(vec1, 20, 10)

  res2 <- which(mat == 1, arr.ind = T)

  #reshuffle
  res <- res[order(res[,1]),]
  res2 <- res2[order(res2[,1]),]

  expect_true(all(res == res2))
})
