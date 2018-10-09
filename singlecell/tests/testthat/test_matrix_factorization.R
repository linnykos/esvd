context("Test matrix factorization")

## .projection_l1 is correct

test_that(".projection_l1 works", {
  set.seed(10)
  current_vec <- rnorm(10)
  other_mat <- matrix(rnorm(50), ncol = 10)
  res <- .projection_l1(current_vec, other_mat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == length(current_vec))
})

test_that(".projection_l1 enforces negativity", {
  set.seed(10)
  current_vec <- rnorm(10)
  other_mat <- matrix(rnorm(50), ncol = 10)
  res <- .projection_l1(current_vec, other_mat)

  expect_true(all(other_mat %*% res <= 0))
})

test_that(".projection_l1 enforces positivity", {
  set.seed(10)
  current_vec <- rnorm(10)
  other_mat <- matrix(rnorm(50), ncol = 10)
  res <- .projection_l1(current_vec, other_mat, direction = ">=")

  expect_true(all(other_mat %*% res >= 0))
})

test_that(".projection_l1 maintains less than 0", {
  trials <- 50

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)

    dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))

    res <- .initialization(dat)
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

    res <- .initialization(dat, family = "gaussian")
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
  res <- .projection_l1(current_vec, other_mat)

  expect_true(sum(abs(res - current_vec)) < 1e-6)
})

test_that(".projection_l1 is actually a projection compared to the all 0 vector", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    current_vec <- abs(rnorm(10))
    other_mat <- -abs(matrix(rnorm(50), ncol = 10))
    new_vec <- .projection_l1(current_vec, other_mat)

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

#########################

## .backtrack_linesearch is correct

test_that(".backtrack_linesearch works", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))

  res <- .initialization(dat)
  u_mat <- res$u_mat
  v_mat <- res$v_mat
  i <- 1
  dat_vec <- dat[i,]
  class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
  grad_vec <- .gradient_vec(dat_vec, u_mat[i,], v_mat)

  res <- .backtrack_linesearch(dat_vec, u_mat[i,], v_mat, grad_vec)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
})

test_that(".backtrack_linesearch actually keeps the negative constraint", {
  trials <- 50

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)

    dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))

    res <- .initialization(dat)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    i <- sample(1:10, 1)

    if(any(!is.na(dat[i,]))){
      dat_vec <- dat[i,]
      class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
      grad_vec <- .gradient_vec(dat_vec, u_mat[i,], v_mat)

      res <- .backtrack_linesearch(dat_vec, u_mat[i,], v_mat, grad_vec)

      u_new <- u_mat[i,] - res*grad_vec
      vec <- v_mat %*% u_new
      all(vec[which(!is.na(dat[i,]))] <= 1e-6)
    } else {
      TRUE
    }
  })

  expect_true(all(bool_vec))
})


test_that(".backtrack_linesearch lowers the objective", {
  trials <- 50

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)

    dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))

    res <- .initialization(dat)
    u_mat <- res$u_mat
    v_mat <- res$v_mat

    i <- sample(1:10, 1)
    if(any(!is.na(dat[i,]))){
      dat_vec <- dat[i,]
      class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
      obj1 <- .evaluate_objective_single(dat_vec, u_mat[i,], v_mat)
      grad_vec <- .gradient_vec(dat_vec, u_mat[i,], v_mat)

      res <- .backtrack_linesearch(dat_vec, u_mat[i,], v_mat, grad_vec)

      u_new <- u_mat[i,] - res*grad_vec
      obj2 <- .evaluate_objective_single(dat_vec, u_new, v_mat)

      obj2 <= obj1 + 1e-6
    } else {
      TRUE
    }
  })

  expect_true(all(bool_vec))
})

##################

## .optimize_row is correct

test_that(".optimize_row works", {
  set.seed(20)
  dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))

  res <- .initialization(dat)
  u_mat <- res$u_mat
  v_mat <- res$v_mat
  i <- 1

  dat_vec <- dat[i,]
  class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
  res <- .optimize_row(dat_vec, u_mat[i,], v_mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == length(u_mat[i,]))
})

test_that(".optimize_row actually lowers the objective", {
  trials <- 25

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))

    res <- .initialization(dat)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    i <- sample(1:10, 1)

    if(any(!is.na(dat[i,]))){
      dat_vec <- dat[i,]
      class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
      u_new <- .optimize_row(dat_vec, u_mat[i,], v_mat)
      obj1 <- .evaluate_objective_single(dat_vec, u_mat[i,], v_mat)
      obj2 <- .evaluate_objective_single(dat_vec, u_new, v_mat)

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

    res <- .initialization(dat)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    j <- sample(1:4, 1)

    if(any(!is.na(dat[,j]))){
      dat_vec <- dat[,j]
      class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
      v_new <- .optimize_row(dat_vec, v_mat[j,], u_mat)
      obj1 <- .evaluate_objective_single(dat_vec, v_mat[j,], u_mat)
      obj2 <- .evaluate_objective_single(dat_vec, v_new, u_mat)

      obj2 <= obj1 + 1e-6
    } else {TRUE}
  })

  expect_true(all(bool_vec))
})

test_that(".optimize_row respects an upper bound", {
  set.seed(20)
  dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))

  res <- .initialization(dat)
  u_mat <- res$u_mat
  v_mat <- res$v_mat
  i <- 1

  dat_vec <- dat[i,]
  class(dat_vec) <- c("exponential", class(dat_vec)[length(class(dat_vec))])
  res1 <- .optimize_row(dat_vec, u_mat[i,], v_mat)
  res2 <- .optimize_row(dat_vec, u_mat[i,], v_mat, max_val = -5)

  expect_true(sum(abs(res1 - res2)) > 1e-6)
  expect_true(all(v_mat %*% res2 >= -5))
})

##################

## .optimize_mat is correct

test_that(".optimize_mat works", {
  set.seed(20)
  dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))
  class(dat) <- c("exponential", class(dat))

  res <- .initialization(dat)
  u_mat <- res$u_mat
  v_mat <- res$v_mat

  res <- .optimize_mat(dat, u_mat, v_mat)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(u_mat)))
})


test_that(".optimize_mat works with parallelization", {
  set.seed(20)
  dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))
  class(dat) <- c("exponential", class(dat))

  res <- .initialization(dat)
  u_mat <- res$u_mat
  v_mat <- res$v_mat

  res1 <- .optimize_mat(dat, u_mat, v_mat, parallelized = F)

  doMC::registerDoMC(cores = 3)
  res2 <- .optimize_mat(dat, u_mat, v_mat, parallelized = T)

  expect_true(sum(abs(res1 - res2)) <= 1e-6)
})

test_that(".optimize_mat keeps the negative constraint", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))
    class(dat) <- c("exponential", class(dat))
    bool <- sample(c(T, F), 1)

    res <- .initialization(dat)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    if(bool){
      res <- .optimize_mat(dat, u_mat, v_mat, bool)
      pred_mat <- res %*% t(v_mat)
    } else {
      res <- .optimize_mat(dat, v_mat, u_mat, bool)
      pred_mat <- u_mat %*% t(res)
    }

    idx <- which(!is.na(dat))

    all(pred_mat[idx] <= -1e-6)
  })

  expect_true(all(bool_vec))
})

test_that(".optimize_mat keeps the positive constraint", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rexp(40), nrow = 10, ncol = 4))
    class(dat) <- c("gaussian", class(dat))
    bool <- sample(c(T, F), 1)

    res <- .initialization(dat, family = "gaussian")
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    if(bool){
      res <- .optimize_mat(dat, u_mat, v_mat, bool)
      pred_mat <- res %*% t(v_mat)
    } else {
      res <- .optimize_mat(dat, v_mat, u_mat, bool)
      pred_mat <- u_mat %*% t(res)
    }

    idx <- which(!is.na(dat))

    all(pred_mat[idx] >= 1e-6)
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

    res <- .initialization(dat)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    obj1 <- .evaluate_objective(dat, u_mat, v_mat)

    if(bool){
      res <- .optimize_mat(dat, u_mat, v_mat, bool)
      obj2 <- .evaluate_objective(dat, res, v_mat)
    } else {
      res <- .optimize_mat(dat, v_mat, u_mat, bool)
      obj2 <- .evaluate_objective(dat, u_mat, res)
    }

    obj2 <= obj1 + 1e-6
  })

  expect_true(all(bool_vec))
})

######################

## .fit_factorization is correct

test_that(".fit_factorization works", {
  set.seed(10)
  dat <- abs(matrix(rexp(20), nrow = 5, ncol = 4))
  init <- .initialization(dat)

  res <- .fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat)

  expect_true(is.list(res))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
  expect_true(ncol(res$u_mat) == 2)
  expect_true(ncol(res$v_mat) == 2)
})

test_that(".fit_factorization works with missing values", {
  set.seed(5)
  dat <- abs(matrix(rexp(20), nrow = 5, ncol = 4))
  init <- .initialization(dat)

  res <- .fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat)

  expect_true(is.list(res))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
  expect_true(ncol(res$u_mat) == 2)
  expect_true(ncol(res$v_mat) == 2)
})

test_that(".fit_factorization preserves the positive entries", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rexp(20), nrow = 5, ncol = 4))
    init <- .initialization(dat)

    res <- .fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat)

    pred_mat <- res$u_mat %*% t(res$v_mat)

    all(pred_mat[which(!is.na(dat))] <= 1e-6)
  })

  expect_true(all(bool_vec))

})

test_that(".fit_factorization can roughly recover the all 1's matrix", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    dat <- abs(matrix(rexp(100), nrow = 10, ncol = 10))
    init <- .initialization(dat)

    fit <- .fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                                          max_iter = 5)

    res1 <- .evaluate_objective(dat, fit$u_mat, fit$v_mat)
    res2 <- .evaluate_objective(dat, matrix(1, ncol = 1, nrow = 10),
                                -matrix(1, ncol = 1, nrow = 10))
    res3 <- .evaluate_objective(dat, abs(matrix(rnorm(10), ncol = 1, nrow = 10)),
                                -abs(matrix(rnorm(10), ncol = 1, nrow = 10)))

    res1 < res2 & res2 < res3
  })
  set.seed(10)

  expect_true(all(bool_vec))

})
