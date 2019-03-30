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

    obj2 < obj1
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

    pred_mat <- abs(matrix(rnorm(100), 10, 10))
    gradient_mat <-  .gradient_mat(dat, pred_mat)

    new_mat <- .adaptive_gradient_step(dat, pred_mat, gradient_mat, k = 2,
                                       direction = direction)

    svd_res <- svd(new_mat)
    all(new_mat > 0) & all(abs(svd_res$d[-(1:2)]) < 1e-3)
  })

  expect_true(all(bool_vec))
})

#####

## .ensure_feasibility is correct

test_that(".ensure_feasibility works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(100), 10, 10))
  direction <- .dictate_direction("gaussian")
  svd_factors <- .svd_projection(dat, k = 2, factors = T)

  res <- .ensure_feasibility(svd_factors$u_mat, svd_factors$v_mat,
                             direction = direction, max_val = 2)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("u_mat", "v_mat")))
})

test_that(".ensure_feasibility respects the constraints", {
  set.seed(10)
  dat <- abs(matrix(rnorm(100), 10, 10))
  direction <- .dictate_direction("gaussian")
  svd_factors <- .svd_projection(dat, k = 2, factors = T)

  res <- .ensure_feasibility(svd_factors$u_mat, svd_factors$v_mat,
                             direction = direction, max_val = 2)

  pred_mat <- res$u_mat %*% t(res$v_mat)

  expect_true(all(pred_mat > 0))
  expect_true(all(pred_mat <= 2))
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
    print(x)
    set.seed(10*x)

    dat <- abs(matrix(rexp(40, 1), nrow = 10, ncol = 4))

    res <- initialization(dat, max_val = -100, max_iter = 5)
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

    obj2 < obj1
  })

  expect_true(all(bool_vec))
})


