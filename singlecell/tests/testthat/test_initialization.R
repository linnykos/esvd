context("Test initialization")

## .matrix_completion is correct

test_that(".matrix_completion works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  res <- .matrix_completion(dat, k = 2)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(dat)))
})

#####

## .svd_projection is correct

test_that(".svd_projection works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  res <- .svd_projection(dat, k = 2)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(dat)))
})

test_that(".svd_projection works with factors", {
  set.seed(10)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  res <- .svd_projection(dat, k = 2, factors = T)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("u_mat", "v_mat")))
  expect_true(all(dim(res$u_mat) == c(nrow(dat), 2)))
  expect_true(all(dim(res$v_mat) == c(ncol(dat), 2)))
})

############

## .adaptive_gradient_step is correct

test_that(".adaptive_gradient_step works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(100), 10, 10))
  class(dat) <- c("gaussian", class(dat)[length(class(dat))])
  direction <- .dictate_direction(class(dat)[1])

  pred_mat <- abs(matrix(rnorm(100), 10, 10))
  gradient_mat <-  .gradient_mat(dat, pred_mat)

  res <- .adaptive_gradient_step(dat, pred_mat, gradient_mat, k = 2,
                                 direction = direction)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(dat)))
})

test_that(".adaptive_gradient_step always decreases the objective", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)

    dat <- abs(matrix(rnorm(100), 10, 10))
    class(dat) <- c("gaussian", class(dat)[length(class(dat))])
    direction <- .dictate_direction(class(dat)[1])

    pred_mat <- abs(matrix(rnorm(100), 10, 10))
    gradient_mat <-  .gradient_mat(dat, pred_mat)

    new_mat <- .adaptive_gradient_step(dat, pred_mat, gradient_mat, k = 2,
                                   direction = direction)

    obj1 <- .evaluate_objective_mat(dat, pred_mat)
    obj2 <- .evaluate_objective_mat(dat, new_mat)

    obj2 <= obj1
  })

  expect_true(all(bool_vec))
})


test_that(".adaptive_gradient_step always gives solutions within the constraint", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)

    dat <- abs(matrix(rnorm(100), 10, 10))
    class(dat) <- c("gaussian", class(dat)[length(class(dat))])
    direction <- .dictate_direction(class(dat)[1])

    pred_mat <- abs(matrix(rnorm(20), 10, 2)) %*% t(abs(matrix(rnorm(20), 10, 2)))
    gradient_mat <-  .gradient_mat(dat, pred_mat)

    new_mat <- .adaptive_gradient_step(dat, pred_mat, gradient_mat, k = 2,
                                       direction = direction)

    svd_res <- svd(new_mat)
    all(new_mat > 0) & all(abs(svd_res$d[-(1:2)]) < 1e-3)
  })

  expect_true(all(bool_vec))
})

############

## .projected_gradient_descent is correct

test_that(".projected_gradient_descent works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(100), 10, 10))
  class(dat) <- c("gaussian", class(dat)[length(class(dat))])
  direction <- .dictate_direction(class(dat)[1])

  res <- .projected_gradient_descent(dat, k = 2, max_iter = 5,
                                     direction = direction)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(dat)))
})

test_that(".projected_gradient_descent decreases the value", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    dat <- abs(matrix(rnorm(100), 10, 10))
    class(dat) <- c("gaussian", class(dat)[length(class(dat))])
    direction <- .dictate_direction(class(dat)[1])

    original_res <- .determine_initial_matrix(dat, class(dat)[1], k = 2)
    res <- .projected_gradient_descent(dat, k = 2, max_iter = 5,
                                       direction = direction)

    obj1 <- .evaluate_objective_mat(dat, original_res)
    obj2 <- .evaluate_objective_mat(dat, res)

    obj2 <= obj1 + 1e-6
  })

  expect_true(all(bool_vec))
})


