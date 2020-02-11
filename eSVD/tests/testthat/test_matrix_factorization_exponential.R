context("Test matrix factorization - Exponential")

# .gradient_vec is correct

test_that(".gradient_vec works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)

  i <- 5
  dat_vec <- dat[i,]
  class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
  res <- .gradient_vec(dat_vec, u_mat[i,], v_mat, n = nrow(dat), p = ncol(dat))

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
  dat_vec <- dat[,j]
  class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
  res <- .gradient_vec(dat_vec, v_mat[j,], u_mat, n = nrow(dat), p = ncol(dat))

  expect_true(is.numeric(res))
  expect_true(length(res) == 2)
})

test_that(".gradient_vec satisfies the gradient definition", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
    u_vec <- abs(rnorm(2))
    u_vec2 <- abs(rnorm(2))
    v_mat <- -abs(matrix(rnorm(8), nrow = 4, ncol = 2))

    i <- sample(1:10, 1)
    dat_vec <- dat[i,]
    class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
    grad <- .gradient_vec(dat_vec, u_vec, v_mat, n = nrow(dat), p = ncol(dat))

    res <- .evaluate_objective_single(dat_vec, u_vec, v_mat, n = nrow(dat), p = ncol(dat))
    res2 <- .evaluate_objective_single(dat_vec, u_vec2, v_mat, n = nrow(dat), p = ncol(dat))

    res2 >= res + as.numeric(grad %*% (u_vec2 - u_vec)) - 1e-6
  })

  expect_true(all(bool_vec))
})

#################


## .evaluate_objective is correct

test_that(".evaluate_objective works", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
  v_mat <- -abs(matrix(rnorm(8), nrow = 4, ncol = 2))
  class(dat) <- c("exponential", class(dat)[length(class(dat))])

  res <- .evaluate_objective(dat, u_mat, v_mat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(!is.nan(res))
})

test_that(".evaluate_objective yields a smaller value under truth", {
  trials <- 100

  avg_obj <- sapply(1:trials, function(x){
    set.seed(x)
    u_mat <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
    v_mat <- -abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    pred_mat <- u_mat %*% t(v_mat)
    dat <- pred_mat
    class(dat) <- c("exponential", class(dat)[length(class(dat))])

    for(i in 1:10){
      for(j in 1:4){
        dat[i,j] <- stats::rexp(1, -pred_mat[i,j])
      }
    }

    res <- .evaluate_objective(dat, u_mat, v_mat)

    u_mat2 <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
    v_mat2 <- -abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    res2 <- .evaluate_objective(dat, u_mat2, v_mat2)

    c(res, res2)
  })

  expect_true(mean(avg_obj[1,]) < mean(avg_obj[2,]))
})

test_that(".evaluate_objective is equal to many .evaluate_objective_single", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
  v_mat <- -abs(matrix(rnorm(8), nrow = 4, ncol = 2))
  class(dat) <- c("exponential", class(dat)[length(class(dat))])

  res <- .evaluate_objective(dat, u_mat, v_mat)

  res2 <- sum(sapply(1:nrow(u_mat), function(x){
    dat_vec <- dat[x,]
    class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
    .evaluate_objective_single(dat[x,], u_mat[x,], v_mat, n = nrow(dat), p = ncol(dat))
  }))

  expect_true(abs(res - res2) <= 1e-6)
})

test_that(".evaluate_objective gives sensible optimal", {
  set.seed(20)
  dat <- abs(matrix(rexp(100, rate = 10), nrow = 10, ncol = 10))
  u_mat <- matrix(2, nrow = 10, ncol = 1)
  v_mat <- -matrix(5, nrow = 10, ncol = 1)
  class(dat) <- c("exponential", class(dat)[length(class(dat))])

  res <- .evaluate_objective(dat, u_mat, v_mat)

  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    u_mat2 <- abs(matrix(rnorm(10, mean = 10), nrow = 10, ncol = 1))
    v_mat2 <- -abs(matrix(rnorm(10, mean = 5), nrow = 10, ncol = 1))
    res2 <- .evaluate_objective(dat, u_mat2, v_mat2)

    res < res2
  })

  expect_true(all(bool_vec))

  u_mat2 <- matrix(50, nrow = 10, ncol = 1)
  v_mat2 <- -matrix(50, nrow = 10, ncol = 1)
  res2 <- .evaluate_objective(dat, u_mat2, v_mat2)
  expect_true(res < res2)

  u_mat2 <- matrix(5, nrow = 10, ncol = 1)
  v_mat2 <- -matrix(5, nrow = 10, ncol = 1)
  res2 <- .evaluate_objective(dat, u_mat2, v_mat2)
  expect_true(res < res2)
})

test_that(".evaluate_objective is correct for rank 1", {
  set.seed(10)
  true_val <- 1/2
  u_mat <- matrix(true_val, nrow = 100, ncol = 1)
  v_mat <- -matrix(true_val, nrow = 100, ncol = 1)
  pred_mat <- u_mat %*% t(v_mat)
  dat <- pred_mat
  class(dat) <- c("exponential", class(dat)[length(class(dat))])

  for(i in 1:nrow(u_mat)){
    for(j in 1:nrow(v_mat)){
      dat[i,j] <- stats::rexp(1, -pred_mat[i,j])
    }
  }

  seq_val <- seq(0.01, 5, length.out = 100)
  nll <- sapply(seq_val, function(x){
    u_mat2 <- matrix(x, nrow = 100, ncol = 1)
    v_mat2 <- -matrix(x, nrow = 100, ncol = 1)
    .evaluate_objective(dat, u_mat2, v_mat2)
  })
  argmin_idx <- which.min(nll)
  min_val <- min(nll)
  supposed_idx <- which.min(abs(seq_val - true_val))
  supposed_val <- nll[supposed_idx]

  expect_true(abs(min_val - supposed_val) <= 1e-6)
})


################

## .evaluate_objective_single is correct

test_that(".evaluate_objective_single works", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
  v_mat <- -abs(matrix(rnorm(8), nrow = 4, ncol = 2))

  i <- 5
  dat_vec <- dat[i,]
  class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
  res <- .evaluate_objective_single(dat_vec, u_mat[i,], v_mat, n = nrow(dat), p = ncol(dat))

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(!is.nan(res))
})

test_that(".evaluate_objective_single yields a smaller value under truth", {
  trials <- 100

  avg_obj <- sapply(1:trials, function(x){
    set.seed(x)
    u_mat <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
    v_mat <- -abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    pred_mat <- u_mat %*% t(v_mat)
    dat <- pred_mat

    for(i in 1:10){
      for(j in 1:4){
        dat[i,j] <- stats::rexp(1, -pred_mat[i,j])
      }
    }

    i <- sample(1:10, 1)
    dat_vec <- dat[i,]
    class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
    res <- .evaluate_objective_single(dat_vec, u_mat[i,], v_mat, n = nrow(dat), p = ncol(dat))

    u_mat2 <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
    v_mat2 <- -abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    res2 <- .evaluate_objective_single(dat_vec, u_mat2[i,], v_mat2, n = nrow(dat), p = ncol(dat))

    c(res, res2)
  })

  expect_true(mean(avg_obj[1,]) < mean(avg_obj[2,]))
})

##########################

test_that("fit_factorization is appropriate for exponential", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    dat <- abs(matrix(rexp(25, 1/2), nrow = 5, ncol = 5))
    class(dat) <- c("exponential", class(dat)[length(class(dat))])
    init <- initialization(dat, family = "exponential")

    fit <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                              max_iter = 5, max_val = -100,
                              family = "exponential")

    res1 <- .evaluate_objective(dat, fit$u_mat, fit$v_mat)
    res2 <- .evaluate_objective(dat, matrix(1, ncol = 1, nrow = 5),
                                -matrix(1, ncol = 1, nrow = 5))
    res3 <- .evaluate_objective(dat, abs(matrix(rnorm(5), ncol = 1, nrow = 5)),
                               -abs(matrix(rnorm(5), ncol = 1, nrow = 5)))

    res1 < res2 & res1 < res3
  })
  set.seed(10)

  expect_true(all(bool_vec))

})

