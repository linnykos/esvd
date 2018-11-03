context("Test matrix factorization - Poisson")

# .gradient_vec is correct

test_that(".gradient_vec works", {
  set.seed(10)
  dat <- abs(matrix(rpois(40, lambda = 1), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- 0.1*matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- 0.1*matrix(rnorm(8), nrow = 4, ncol = 2)

  i <- 5
  dat_vec <- dat[i,]
  class(dat_vec) <- c("poisson", class(dat_vec)[length(class(dat_vec))])
  res <- .gradient_vec(dat_vec, u_mat[i,], v_mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == 2)
})

test_that(".gradient_vec works for the other direction", {
  set.seed(8)
  dat <- abs(matrix(rpois(40, lambda = 1), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- 0.1*matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- 0.1*matrix(rnorm(8), nrow = 4, ncol = 2)

  j <- 2
  dat_vec <- dat[,j]
  class(dat_vec) <- c("poisson", class(dat_vec)[length(class(dat_vec))])
  res <- .gradient_vec(dat_vec, v_mat[j,], u_mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == 2)
})

test_that(".gradient_vec satisfies the gradient definition", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    dat <- abs(matrix(rpois(40, lambda = 1), nrow = 10, ncol = 4))
    u_vec <- 0.1*abs(rnorm(2))
    u_vec2 <- 0.1*abs(rnorm(2))
    v_mat <- 0.1*abs(matrix(rnorm(8), nrow = 4, ncol = 2))

    i <- sample(1:10, 1)
    dat_vec <- dat[i,]
    class(dat_vec) <- c("poisson", class(dat_vec)[length(class(dat_vec))])
    grad <- .gradient_vec(dat_vec, u_vec, v_mat)

    res <- .evaluate_objective_single(dat_vec, u_vec, v_mat)
    res2 <- .evaluate_objective_single(dat_vec, u_vec2, v_mat)

    res2 >= res + as.numeric(grad %*% (u_vec2 - u_vec)) - 1e-6
  })

  expect_true(all(bool_vec))
})

#################


## .evaluate_objective is correct

test_that(".evaluate_objective works", {
  set.seed(20)
  dat <- abs(matrix(rpois(40, lambda = 0.1), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- 0.1*abs(matrix(rnorm(20), nrow = 10, ncol = 2))
  v_mat <- 0.1*abs(matrix(rnorm(8), nrow = 4, ncol = 2))
  class(dat) <- c("poisson", class(dat)[length(class(dat))])

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
    u_mat <- 0.1*abs(matrix(rnorm(20), nrow = 10, ncol = 2))
    v_mat <- 0.1*abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    pred_mat <- u_mat %*% t(v_mat)
    dat <- pred_mat
    class(dat) <- c("poisson", class(dat)[length(class(dat))])

    for(i in 1:10){
      for(j in 1:4){
        dat[i,j] <- stats::rpois(1, lambda = exp(pred_mat[i,j]))
      }
    }

    res <- .evaluate_objective(dat, u_mat, v_mat)

    u_mat2 <- 0.1*abs(matrix(rnorm(20), nrow = 10, ncol = 2))
    v_mat2 <- 0.1*abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    res2 <- .evaluate_objective(dat, u_mat2, v_mat2)

    c(res, res2)
  })

  expect_true(mean(avg_obj[1,]) < mean(avg_obj[2,]))
})

test_that(".evaluate_objective is equal to many .evaluate_objective_single", {
  set.seed(20)
  dat <- matrix(rpois(40, lambda = 1), nrow = 10, ncol = 4)
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- 0.1*abs(matrix(rnorm(20), nrow = 10, ncol = 2))
  v_mat <- 0.1*abs(matrix(rnorm(8), nrow = 4, ncol = 2))
  class(dat) <- c("poisson", class(dat)[length(class(dat))])

  res <- .evaluate_objective(dat, u_mat, v_mat)

  res2 <- sum(sapply(1:nrow(u_mat), function(x){
    dat_vec <- dat[x,]
    class(dat_vec) <- c("poisson", class(dat_vec)[length(class(dat_vec))])
    .evaluate_objective_single(dat_vec, u_mat[x,], v_mat)
  }))

  expect_true(abs(res - res2) <= 1e-6)
})

test_that(".evaluate_objective gives sensible optimal", {
  set.seed(20)
  dat <- abs(matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10))
  u_mat <- matrix(sqrt(log(5)), nrow = 10, ncol = 1)
  v_mat <- matrix(sqrt(log(5)), nrow = 10, ncol = 1)
  class(dat) <- c("poisson", class(dat)[length(class(dat))])

  res <- .evaluate_objective(dat, u_mat, v_mat)

  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    u_mat2 <- 0.5*abs(matrix(rnorm(10, mean = 10), nrow = 10, ncol = 1))
    v_mat2 <- 0.5*abs(matrix(rnorm(10, mean = 5), nrow = 10, ncol = 1))
    res2 <- .evaluate_objective(dat, u_mat2, v_mat2)

    res < res2
  })

  expect_true(all(bool_vec))

  u_mat2 <- matrix(2, nrow = 10, ncol = 1)
  v_mat2 <- matrix(2, nrow = 10, ncol = 1)
  res2 <- .evaluate_objective(dat, u_mat2, v_mat2)
  expect_true(res < res2)
})


################

## .evaluate_objective_single is correct

test_that(".evaluate_objective_single works", {
  set.seed(20)
  dat <- abs(matrix(rpois(40, lambda = 1), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- 0.1*abs(matrix(rnorm(20), nrow = 10, ncol = 2))
  v_mat <- 0.1*abs(matrix(rnorm(8), nrow = 4, ncol = 2))

  i <- 5
  dat_vec <- dat[i,]
  class(dat_vec) <- c("poisson", class(dat_vec)[length(class(dat_vec))])
  res <- .evaluate_objective_single(dat_vec, u_mat[i,], v_mat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(!is.nan(res))
})

test_that(".evaluate_objective_single yields a smaller value under truth", {
  trials <- 100

  avg_obj <- sapply(1:trials, function(x){
    set.seed(x)
    u_mat <- 0.1*abs(matrix(rnorm(20), nrow = 10, ncol = 2))
    v_mat <- 0.1*abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    pred_mat <- u_mat %*% t(v_mat)
    dat <- pred_mat

    for(i in 1:10){
      for(j in 1:4){
        dat[i,j] <- stats::rpois(1, lambda = exp(pred_mat[i,j]))
      }
    }

    i <- sample(1:10, 1)
    dat_vec <- dat[i,]
    class(dat_vec) <- c("poisson", class(dat_vec)[length(class(dat_vec))])
    res <- .evaluate_objective_single(dat_vec, u_mat[i,], v_mat)

    u_mat2 <- 0.1*abs(matrix(rnorm(20), nrow = 10, ncol = 2))
    v_mat2 <- 0.1*abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    res2 <- .evaluate_objective_single(dat_vec, u_mat2[i,], v_mat2)

    c(res, res2)
  })

  expect_true(mean(avg_obj[1,]) < mean(avg_obj[2,]))
})
