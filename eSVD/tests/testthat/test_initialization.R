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

test_that(".svd_projection correctly returns the matrix if the matrix is already low-rank", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    u_mat <- matrix(rnorm(50), 25, 2)
    v_mat <- matrix(rnorm(60), 30, 2)
    dat <- u_mat %*% t(v_mat)

    res <- .svd_projection(dat, k = 2, factors = F)

    sum(abs(res - dat)) <= 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".svd_projection works on a difficult example", {
  load("../assets/svd1.RData")
  res <- .svd_projection(mat, k = 2)

  expect_true(all(dim(mat) == dim(res)))
})

test_that(".svd_projection works on another difficult example", {
  load("../assets/svd2.RData")
  expect_true(Matrix::rankMatrix(mat) == 2)
  res <- .svd_projection(mat, k = 2)

  expect_true(sum(abs(mat - res)) <= 1e-6)
})

############

## .adaptive_gradient_step is correct

test_that(".adaptive_gradient_step works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(100), 10, 10))
  class(dat) <- c("gaussian", class(dat)[length(class(dat))])
  direction <- .dictate_direction(class(dat)[1])

  pred_mat <- -abs(matrix(rnorm(100), 10, 10))
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

    pred_mat <- -abs(matrix(rnorm(100), 10, 10))
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

    pred_mat <- -abs(matrix(rnorm(20), 10, 2)) %*% t(abs(matrix(rnorm(20), 10, 2)))
    gradient_mat <-  .gradient_mat(dat, pred_mat)

    new_mat <- .adaptive_gradient_step(dat, pred_mat, gradient_mat, k = 2,
                                       direction = direction)

    svd_res <- svd(new_mat)
    all(new_mat < 0) & all(abs(svd_res$d[-(1:2)]) < 1e-3)
  })

  expect_true(all(bool_vec))
})

############

## .projected_gradient_descent is correct

test_that(".projected_gradient_descent works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(100), 10, 10))
  class(dat) <- c("curved_gaussian", class(dat)[length(class(dat))])
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
    class(dat) <- c("curved_gaussian", class(dat)[length(class(dat))])
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

#####################

## initialization is correct

test_that("initialization works", {
  set.seed(10)
  true_val <- 1/2
  u_mat <- matrix(true_val, nrow = 100, ncol = 1)
  v_mat <- matrix(true_val, nrow = 100, ncol = 1)
  pred_mat <- u_mat %*% t(v_mat)
  dat <- pred_mat
  class(dat) <- c("curved_gaussian", class(dat)[length(class(dat))])

  for(i in 1:nrow(u_mat)){
    for(j in 1:nrow(v_mat)){
      dat[i,j] <- abs(stats::rnorm(1, mean = 1/pred_mat[i,j], sd = 1/(2*pred_mat[i,j])))
    }
  }

  res <- initialization(dat, k = 2, family = "gaussian", max_val = 100)

  expect_true(is.list(res))
})

test_that("initialization respects max_val", {
  trials <- 25

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- matrix(rexp(40), nrow = 10, ncol = 4)

    res <- initialization(dat, max_val = 100, family = "exponential")
    u_mat <- res$u_mat
    v_mat <- res$v_mat

    all(abs(u_mat %*% t(v_mat)) <= 100)
  })

  expect_true(all(bool_vec))
})


test_that("initialization is meaningful for a tricky instance", {
  load("../assets/initialization1.RData")

  res <- initialization(dat, family = "curved_gaussian", k = 2, max_val = 100,
                        scalar = 4)

  pred_mat <- res$u_mat %*% t(res$v_mat)
  expect_true(length(unique(pred_mat)) >= prod(dim(pred_mat))/2)
})

#########################

## .project_rank_feasibility is correct

test_that(".project_rank_feasibility works", {
  set.seed(10)
  cov_mat <- matrix(0, 10, 10)
  cov_mat[1:5,1:5] <- 0.3
  cov_mat[6:10,6:10] <- 0.5
  diag(cov_mat) <- 1

  dat <- MASS::mvrnorm(n = 50, mu = rep(2, 10), Sigma = cov_mat)

  res <- .project_rank_feasibility(dat, k = 2, direction = ">=")

  expect_true(all(dim(res$matrix) == dim(dat)))
  expect_true(Matrix::rankMatrix(res$matrix) == 2)
  expect_true(all(res$matrix >= 0))
})

test_that(".project_rank_feasibility is able to iterate more than twice", {
  trials <- 200

  bool_vec <- sapply(1:trials, function(i){
    set.seed(i)
    obs <- matrix(abs(rnorm(25)), 5, 5)
    res <- .project_rank_feasibility(obs, k = 2, direction = "<=")
    res$iter >= 2
  })

  expect_true(any(bool_vec))
})

test_that(".project_rank_feasibility can be the same as SVD in simple settings",{
  set.seed(10)
  dat <- matrix(abs(rnorm(25, mean = 10)), 5, 5)
  res <-.project_rank_feasibility(dat, k = 2, direction = ">=")

  expect_true(res$iter == 1)
})

test_that(".project_rank_feasibility is doing something meaningful", {
  # generate a bunch of matrices and their solutions, and see if the approximated
  #   matrices are indeed the intended solutions
  trials <- 100

  dat_list <- lapply(1:trials, function(i){
    set.seed(i)
    matrix(abs(rnorm(25, mean = 1)), 5, 5) # we use a high signal since our method doesn't seem to do that well otherwise...
  })

  res_list <- lapply(dat_list, function(x){
    res <- .project_rank_feasibility(x, k = 2, direction = ">=")
    res
  })

  iter_vec <- sapply(res_list, function(x){x$iter})
  expect_true(any(iter_vec >= 2))
  res_list <- lapply(res_list, function(x){x$matrix})

  # compute the cross-error
  error_mat <- sapply(1:trials, function(i){
    vec <- sapply(1:trials, function(j){
      .l2norm(dat_list[[i]] - res_list[[j]])
    })
    vec/min(vec)
  })

  expect_true(all(diag(error_mat) <= 1.01))
})

test_that(".project_rank_feasibility handles non-convergence settings gracefully (special example)", {
  set.seed(1)
  n <- 50
  u_mat <- matrix(abs(rnorm(50)), ncol = 1)
  v_mat <- matrix(abs(rnorm(50)), ncol = 1)
  pred_mat <- u_mat %*% t(v_mat)

  dat <- pred_mat

  for(i in 1:10){
    for(j in 1:4){
      dat[i,j] <- abs(stats::rnorm(1, mean = 1/pred_mat[i,j], sd = 1/(2*pred_mat[i,j])))
    }
  }

  family = "curved_gaussian"
  k = 2
  tol = 1e-3
  max_val = 100
  max_iter = 10
  verbose = F

  direction <- .dictate_direction(family)
  if(!is.na(direction) & !is.na(max_val)){
    stopifnot((direction == ">=" & max_val > 0) | (direction == "<=" & max_val < 0))
  }

  # initialize
  dat <- .matrix_completion(dat, k = k)
  if(length(class(dat)) == 1) class(dat) <- c(family, class(dat)[length(class(dat))])

  min_val <- min(dat[which(dat > 0)])
  dat[which(dat <= 0)] <- min_val/2
  pred_mat <- .mean_transformation(dat, family)
  direction <- .dictate_direction(family)

  res <- .project_rank_feasibility(pred_mat, k = k, direction = direction, max_val = max_val,
                                   max_iter = max_iter)

  expect_true(Matrix::rankMatrix(res$matrix) == 2)
  expect_true(all(res$matrix >= 0))
  expect_true(all(res$matrix <= max_val))
})

###########################3

## .absolute_threshold is correct

test_that(".absolute_threshold works", {
  mat <- matrix(1:25, 5, 5)
  res <- .absolute_threshold(mat, direction = ">=", max_val = 10)

  expect_true(all(dim(mat) == dim(res)))
  expect_true(all(res <= 10))
})

test_that(".absolute_threshold correctly sets things below the threshold", {
  mat <- -matrix(abs(rnorm(100)), 10, 10)
  res <- .absolute_threshold(mat, direction = "<=", max_val = -2.5)

  expect_true(all(res >= -2.5))
})
