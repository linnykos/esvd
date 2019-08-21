context("Test matrix factorization")
library(NMF)

## .optimize_row is correct

test_that(".optimize_row works", {
  set.seed(20)
  dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))

  res <- initialization(dat)
  u_mat <- res$u_mat
  v_mat <- res$v_mat
  i <- 1

  dat_vec <- dat[i,]
  class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
  res <- .optimize_row(dat_vec, u_mat[i,], v_mat, max_val = -100, n = nrow(dat), p = ncol(dat))

  expect_true(is.numeric(res))
  expect_true(length(res) == length(u_mat[i,]))
})

test_that(".optimize_row actually lowers the objective", {
  trials <- 25

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))

    res <- initialization(dat, max_val = -100)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    i <- sample(1:10, 1)

    if(any(!is.na(dat[i,]))){
      dat_vec <- dat[i,]
      class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
      u_new <- .optimize_row(dat_vec, u_mat[i,], v_mat, max_val = -100, n = nrow(dat), p = ncol(dat))
      obj1 <- .evaluate_objective_single(dat_vec, u_mat[i,], v_mat, n = nrow(dat), p = ncol(dat))
      obj2 <- .evaluate_objective_single(dat_vec, u_new, v_mat, n = nrow(dat), p = ncol(dat))

      obj2 <= obj1 + 1e-6
    } else {TRUE}
  })

  expect_true(all(bool_vec))
})

test_that(".optimize_row works the other way", {
  trials <- 25

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))

    res <- initialization(dat, max_val = -100)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    j <- sample(1:4, 1)

    if(any(!is.na(dat[,j]))){
      dat_vec <- dat[,j]
      class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
      v_new <- .optimize_row(dat_vec, v_mat[j,], u_mat, max_val = -100, n = nrow(dat), p = ncol(dat))
      obj1 <- .evaluate_objective_single(dat_vec, v_mat[j,], u_mat, n = nrow(dat), p = ncol(dat))
      obj2 <- .evaluate_objective_single(dat_vec, v_new, u_mat, n = nrow(dat), p = ncol(dat))

      obj2 <= obj1 + 1e-6
    } else {TRUE}
  })

  expect_true(all(bool_vec))
})

test_that(".optimize_row respects an upper bound", {
  set.seed(20)
  dat <- matrix(rexp(40), nrow = 10, ncol = 4)

  res <- initialization(dat, max_val = -5, family = "exponential")
  u_mat <- res$u_mat
  v_mat <- res$v_mat
  i <- 1

  dat_vec <- dat[i,]
  class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
  res1 <- .optimize_row(dat_vec, u_mat[i,], v_mat, max_val = -100, n = nrow(dat), p = ncol(dat))
  res2 <- .optimize_row(dat_vec, u_mat[i,], v_mat, max_val = -5, n = nrow(dat), p = ncol(dat))

  expect_true(sum(abs(res1 - res2)) > 1e-6)
  expect_true(all(v_mat %*% res2 >= -5))
})

##################

## .optimize_mat is correct

test_that(".optimize_mat works", {
  set.seed(20)
  dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))
  class(dat) <- c("exponential", class(dat))

  res <- initialization(dat, max_val = -100)
  u_mat <- res$u_mat
  v_mat <- res$v_mat

  res <- .optimize_mat(dat, u_mat, v_mat, max_val = -100)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(u_mat)))
})


test_that(".optimize_mat works with parallelization", {
  set.seed(20)
  dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))
  class(dat) <- c("exponential", class(dat))

  res <- initialization(dat, max_val = -100)
  u_mat <- res$u_mat
  v_mat <- res$v_mat

  res1 <- .optimize_mat(dat, u_mat, v_mat, parallelized = F, max_val = -100)

  doMC::registerDoMC(cores = 2)
  res2 <- .optimize_mat(dat, u_mat, v_mat, parallelized = T, max_val = -100)

  expect_true(sum(abs(res1 - res2)) <= 1e-6)
})

test_that(".optimize_mat keeps the negative constraint", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))
    class(dat) <- c("exponential", class(dat))
    bool <- sample(c(T, F), 1)

    res <- initialization(dat, max_val = -100)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    if(bool){
      res <- .optimize_mat(dat, u_mat, v_mat, bool, max_val = -100)
      pred_mat <- res %*% t(v_mat)
    } else {
      res <- .optimize_mat(dat, v_mat, u_mat, bool, max_val = -100)
      pred_mat <- u_mat %*% t(res)
    }

    idx <- which(!is.na(dat))

    all(pred_mat[idx] <= -1e-6)
  })

  expect_true(all(bool_vec))
})

test_that(".optimize_mat keeps the negative constraint for Gaussians", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))
    class(dat) <- c("gaussian", class(dat))
    bool <- sample(c(T, F), 1)

    res <- initialization(dat, family = "gaussian", max_val = -100)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    if(bool){
      res <- .optimize_mat(dat, u_mat, v_mat, bool, max_val = -100)
      pred_mat <- res %*% t(v_mat)
    } else {
      res <- .optimize_mat(dat, v_mat, u_mat, bool, max_val = -100)
      pred_mat <- u_mat %*% t(res)
    }

    idx <- which(!is.na(dat))

    all(pred_mat[idx] <= -1e-6)
  })

  expect_true(all(bool_vec))
})


test_that(".optimize_mat lowers the objective value", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))
    class(dat) <- c("exponential", class(dat))
    bool <- sample(c(T, F), 1)

    res <- initialization(dat, max_val = -100)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    obj1 <- .evaluate_objective(dat, u_mat, v_mat)

    if(bool){
      res <- .optimize_mat(dat, u_mat, v_mat, bool, max_val = -100)
      obj2 <- .evaluate_objective(dat, res, v_mat)
    } else {
      res <- .optimize_mat(dat, v_mat, u_mat, bool, max_val = -100)
      obj2 <- .evaluate_objective(dat, u_mat, res)
    }

    obj2 <= obj1 + 1e-6
  })

  expect_true(all(bool_vec))
})

#######################

## .frank_wolfe is correct

test_that(".frank_wolfe is able to solve the following LP", {
  load("../assets/frank_wolfe1.RData")
  res <- .frank_wolfe(grad_vec, other_mat, other_bound = -150)

  expect_true(is.numeric(res))
  expect_true(length(res) == 2)
  expect_true(all(other_mat %*% res >= -150-1e-3))
  expect_true(all(other_mat %*% res <= 0))
})

######################

## fit_factorization is correct

test_that("fit_factorization works", {
  set.seed(10)
  dat <- abs(matrix(rexp(20), nrow = 5, ncol = 4))
  init <- initialization(dat, max_val = -100)

  res <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                            max_val = -100)

  expect_true(is.list(res))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
  expect_true(ncol(res$u_mat) == 2)
  expect_true(ncol(res$v_mat) == 2)
})

test_that("fit_factorization works for Gaussian with scalar setting", {
  set.seed(10)
  dat <- abs(matrix(rexp(20), nrow = 5, ncol = 4))

  init <- initialization(dat, family = "gaussian", max_val = -100)
  res <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                            family = "gaussian",
                            max_val = -100, scalar = 3)

  expect_true(is.list(res))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
  expect_true(ncol(res$u_mat) == 2)
  expect_true(ncol(res$v_mat) == 2)
  expect_true(min(res$u_mat %*% t(res$v_mat)) >= -100)
  expect_true(max(res$u_mat %*% t(res$v_mat)) <= 0)
})

test_that("fit_factorization for Gaussian with scalar setting respects max_val", {
  set.seed(10)
  dat <- abs(matrix(rexp(20), nrow = 5, ncol = 4))

  init <- initialization(dat, family = "gaussian", max_val = -5)
  res <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                            family = "gaussian",
                            max_val = -5, scalar = 100)

  expect_true(is.list(res))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
  expect_true(ncol(res$u_mat) == 2)
  expect_true(ncol(res$v_mat) == 2)
  expect_true(min(res$u_mat %*% t(res$v_mat)) >= -5)
  expect_true(max(res$u_mat %*% t(res$v_mat)) <= 0)
})

test_that("fit_factorization works with missing values", {
  set.seed(5)
  dat <- abs(matrix(rexp(20), nrow = 5, ncol = 4))
  dat[sample(prod(dim(dat)), 5)] <- NA
  init <- initialization(dat, max_val = -100)

  res <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                            max_val = -100)

  expect_true(is.list(res))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
  expect_true(ncol(res$u_mat) == 2)
  expect_true(ncol(res$v_mat) == 2)
})

test_that("fit_factorization preserves the negative entries", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rexp(20), nrow = 5, ncol = 4))
    init <- initialization(dat, max_val = -100)

    res <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                              max_val = -100)

    pred_mat <- res$u_mat %*% t(res$v_mat)

    all(pred_mat[which(!is.na(dat))] <= 1e-6)
  })

  expect_true(all(bool_vec))

})

test_that("fit_factorization works with Gaussian and missing values", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rnorm(20), nrow = 5, ncol = 4))
    dat[sample(1:prod(dim(dat)),3)] <- NA
    init <- initialization(dat, max_val = -100, family = "gaussian")

    res <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                              max_val = -100, max_iter = 3, family = "gaussian")

    pred_mat <- res$u_mat %*% t(res$v_mat)

    all(pred_mat >= -100+1e-6) & all(pred_mat <= 0)
  })

  expect_true(all(bool_vec))
})

test_that("fit_factorization gives similar results if only one value is missing", {
  set.seed(10)
  dat <- abs(matrix(rexp(100), nrow = 10, ncol = 10))

  dat2 <- dat; dat2[2,1] <- NA
  dat3 <- dat
  for(i in 1:nrow(dat3)){
    dat3[i, sample(1:ncol(dat3), 2)] <- NA
  }

  init <- initialization(dat, max_val = -100, family = "gaussian")
  res <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                            max_val = -100, max_iter = 100, family = "gaussian")

  init2 <- initialization(dat2, max_val = -100, family = "gaussian")
  res2 <- fit_factorization(dat2, u_mat = init2$u_mat, v_mat = init2$v_mat,
                            max_val = -100, max_iter = 100, family = "gaussian")

  init3 <- initialization(dat3, max_val = -100, family = "gaussian")
  res3 <- fit_factorization(dat3, u_mat = init3$u_mat, v_mat = init3$v_mat,
                             max_val = -100, max_iter = 100, family = "gaussian")

  expect_true(.l2norm(res$u_mat - res2$u_mat) < .l2norm(res$u_mat - res3$u_mat))
})

test_that("fit_factorization respects contraints for all values, even the missing ones", {
  set.seed(10)
  dat <- abs(matrix(rexp(100), nrow = 10, ncol = 10))
  for(i in 1:nrow(dat)){
    dat[i, sample(1:ncol(dat), 2)] <- NA
  }

  init <- initialization(dat, max_val = -100, family = "gaussian")
  res <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                             max_val = -100, max_iter = 25, family = "gaussian")

  pred_mat <- res$u_mat %*% t(res$v_mat)

  expect_true(all(pred_mat < 0))
})


test_that("fit_factorization can roughly recover the all 1's matrix", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    dat <- abs(matrix(rexp(100), nrow = 10, ncol = 10))
    init <- initialization(dat, max_val = -100)

    fit <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                                          max_iter = 5, max_val = -100)

    res1 <- .evaluate_objective(dat, fit$u_mat, fit$v_mat)
    res2 <- .evaluate_objective(dat, matrix(1, ncol = 1, nrow = 10),
                                -matrix(1, ncol = 1, nrow = 10))
    res3 <- .evaluate_objective(dat, abs(matrix(rnorm(10), ncol = 1, nrow = 10)),
                                -abs(matrix(rnorm(10), ncol = 1, nrow = 10)))

    res1 < res2 & res1 < res3
  })
  set.seed(10)

  expect_true(all(bool_vec))
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

test_that(".reparameterize preserves the inner products", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    u_mat <- MASS::mvrnorm(5, rep(0, 5), diag(5))
    v_mat <- MASS::mvrnorm(5, rep(1, 5), 2*diag(5))

    res <- .reparameterize(u_mat, v_mat)

    pred_mat <- u_mat %*% t(v_mat)
    pred_mat2 <- res$u_mat %*% t(res$v_mat)

    sum(abs(pred_mat - pred_mat2)) <= 1e-6
  })

  expect_true(all(bool_vec))
})

