context("Test matrix factorization - Gaussian")

# .gradient_vec is correct

test_that(".gradient_vec works", {
  set.seed(10)
  dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)

  i <- 5
  dat_vec <- dat[i,]
  attr(dat_vec, "family") <- "gaussian"
  res <- .gradient_vec(dat_vec, u_mat[i,], v_mat, n = nrow(dat), p = ncol(dat))

  expect_true(is.numeric(res))
  expect_true(length(res) == 2)
})


test_that(".gradient_vec works for the other direction", {
  set.seed(8)
  dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)

  j <- 2
  dat_vec <- dat[,j]
  attr(dat_vec, "family") <- "gaussian"
  res <- .gradient_vec(dat_vec, v_mat[j,], u_mat, n = nrow(dat), p = ncol(dat))

  expect_true(is.numeric(res))
  expect_true(length(res) == 2)
})

test_that(".gradient_vec satisfies the gradient definition", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
    u_vec <- rnorm(2)
    u_vec2 <- rnorm(2)
    v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)

    i <- sample(1:10, 1)
    dat_vec <- dat[i,]
    attr(dat_vec, "family") <- "gaussian"
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
  dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)
  attr(dat, "family") <- "gaussian"

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
    u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
    v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)
    pred_mat <- u_mat %*% t(v_mat)
    dat <- pred_mat
    attr(dat, "family") <- "gaussian"

    for(i in 1:10){
      for(j in 1:4){
        dat[i,j] <- stats::rnorm(1, mean = pred_mat[i,j])
      }
    }

    res <- .evaluate_objective(dat, u_mat, v_mat)

    u_mat2 <- matrix(rnorm(20), nrow = 10, ncol = 2)
    v_mat2 <- matrix(rnorm(8), nrow = 4, ncol = 2)
    res2 <- .evaluate_objective(dat, u_mat2, v_mat2)

    c(res, res2)
  })

  expect_true(mean(avg_obj[1,]) < mean(avg_obj[2,]))
})

test_that(".evaluate_objective is correct for rank 1", {
  set.seed(10)
  true_val <- 1/2
  u_mat <- matrix(true_val, nrow = 100, ncol = 1)
  v_mat <- matrix(true_val, nrow = 100, ncol = 1)
  pred_mat <- u_mat %*% t(v_mat)
  dat <- pred_mat
  attr(dat, "family") <- "gaussian"

  for(i in 1:nrow(u_mat)){
    for(j in 1:nrow(v_mat)){
      dat[i,j] <- stats::rnorm(1, mean = pred_mat[i,j])
    }
  }

  seq_val <- seq(0.01, 5, length.out = 100)
  nll <- sapply(seq_val, function(x){
    u_mat2 <- matrix(x, nrow = 100, ncol = 1)
    v_mat2 <- matrix(x, nrow = 100, ncol = 1)
    .evaluate_objective(dat, u_mat2, v_mat2)
  })

  min_val <- seq_val[which.min(nll)]
  supposed_val <- seq_val[which.min(abs(seq_val - true_val))]

  expect_true(abs(min_val - supposed_val) <= 1e-6)
})


test_that(".evaluate_objective is equal to many .evaluate_objective_single", {
  set.seed(20)
  dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)
  attr(dat, "family") <- "gaussian"

  res <- .evaluate_objective(dat, u_mat, v_mat)

  res2 <- sum(sapply(1:nrow(u_mat), function(x){
    dat_vec <- dat[x,]
    attr(dat_vec, "family") <- "gaussian"
    .evaluate_objective_single(dat_vec, u_mat[x,], v_mat, n = nrow(dat), p = ncol(dat))
  }))

  expect_true(abs(res - res2) <= 1e-6)
})

test_that(".evaluate_objective gives sensible optimal", {
  set.seed(20)
  dat <- matrix(rnorm(100, mean = 1/10), nrow = 10, ncol = 10)
  u_mat <- matrix(1/2, nrow = 10, ncol = 1)
  v_mat <- matrix(1/5, nrow = 10, ncol = 1)
  attr(dat, "family") <- "gaussian"

  res <- .evaluate_objective(dat, u_mat, v_mat)

  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    u_mat2 <- matrix(rnorm(10, mean = 10), nrow = 10, ncol = 1)
    v_mat2 <- matrix(rnorm(10, mean = 5), nrow = 10, ncol = 1)
    res2 <- .evaluate_objective(dat, u_mat2, v_mat2)

    res < res2
  })

  expect_true(all(bool_vec))

  u_mat2 <- matrix(50, nrow = 10, ncol = 1)
  v_mat2 <- matrix(50, nrow = 10, ncol = 1)
  res2 <- .evaluate_objective(dat, u_mat2, v_mat2)
  expect_true(res < res2)

  u_mat2 <- matrix(5, nrow = 10, ncol = 1)
  v_mat2 <- matrix(5, nrow = 10, ncol = 1)
  res2 <- .evaluate_objective(dat, u_mat2, v_mat2)
  expect_true(res < res2)
})


################

## .evaluate_objective_single is correct

test_that(".evaluate_objective_single works", {
  set.seed(20)
  dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)

  i <- 5
  dat_vec <- dat[i,]
  attr(dat_vec, "family") <- "gaussian"
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
    u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
    v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)
    pred_mat <- u_mat %*% t(v_mat)
    dat <- pred_mat

    for(i in 1:10){
      for(j in 1:4){
        dat[i,j] <- stats::rnorm(1, pred_mat[i,j])
      }
    }

    i <- sample(1:10, 1)
    dat_vec <- dat[i,]
    attr(dat_vec, "family") <- "gaussian"
    res <- .evaluate_objective_single(dat_vec, u_mat[i,], v_mat, n = nrow(dat), p = ncol(dat))

    u_mat2 <- matrix(rnorm(20), nrow = 10, ncol = 2)
    v_mat2 <- matrix(rnorm(8), nrow = 4, ncol = 2)
    res2 <- .evaluate_objective_single(dat_vec, u_mat2[i,], v_mat2, n = nrow(dat), p = ncol(dat))

    c(res, res2)
  })

  expect_true(mean(avg_obj[1,]) < mean(avg_obj[2,]))
})

#########

## .evaluate_objective_mat.gaussian is correct

test_that(".evaluate_objective_mat.gaussian works", {
  set.seed(20)
  dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)
  pred_mat <- u_mat %*% t(v_mat)
  attr(dat, "family") <- "gaussian"

  res <- .evaluate_objective_mat(dat, pred_mat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(!is.nan(res))
})

test_that(".evaluate_objective_mat is the same as .evaluate_objective", {
  set.seed(10)
  dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)
  pred_mat <- u_mat %*% t(v_mat)
  attr(dat, "family") <- "gaussian"

  res <- .evaluate_objective_mat(dat, pred_mat)
  res2 <- .evaluate_objective(dat, u_mat, v_mat)

  expect_true(abs(res - res2) <= 1e-6)
})


#########

## .gradient_mat.gaussian is correct

test_that(".gradient_mat.gaussian works", {
  set.seed(20)
  dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
  u_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  v_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)
  pred_mat <- u_mat %*% t(v_mat)
  attr(dat, "family") <- "gaussian"

  res <- .gradient_mat(dat, pred_mat)

  expect_true(is.numeric(res))
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(dat)))
})

test_that(".gradient_mat.gaussiann is a proper gradient", {
  trials <- 500

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    dat <- matrix(rnorm(40), nrow = 10, ncol = 4)
    pred_mat <- matrix(rnorm(40), nrow = 10, ncol = 4)
    pred_mat2 <- matrix(rnorm(40), nrow = 10, ncol = 4)

    attr(dat, "family") <- "gaussian"
    grad <-  .gradient_mat(dat, pred_mat)

    res <- .evaluate_objective_mat(dat, pred_mat)
    res2 <- .evaluate_objective_mat(dat, pred_mat2)

    res2 >= res + t(as.numeric(grad))%*%as.numeric(pred_mat2 - pred_mat) - 1e-6
  })

  expect_true(all(bool_vec))
})


#######################

### Gaussian specific tests: initialization


test_that("fit_factorization is about the (roughly) same as SVD", {
  trials <- 10

  diff_val <- sapply(1:trials, function(x){
    set.seed(13*x)
    dat <- matrix(abs(rnorm(25, 2, 1)), nrow = 5, ncol = 5)
    attr(dat, "family") <- "gaussian"
    init <- initialization(dat, family = "gaussian", max_val = 100)

    pred_mat1 <- init$u_mat %*% t(init$v_mat)

    svd_res <- svd(dat)
    pred_mat2 <- svd_res$u[,1:2]%*%diag(svd_res$d[1:2])%*%t(svd_res$v[,1:2])

    sum(abs(pred_mat1 - pred_mat2))
  })

  expect_true(sum(abs(diff_val)) <= 1e-6)
})




###########################3

### Gaussian specific tests: fit_factorzation


test_that("fit_factorization is appropriate for gaussians", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    dat <- matrix(abs(rnorm(25, 2, 1)), nrow = 5, ncol = 5)
    attr(dat, "family") <- "gaussian"
    init <- initialization(dat, family = "gaussian", max_val = 100)

    fit <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                             max_iter = 5, max_val = 100,
                             family = "gaussian")

    res1 <- .evaluate_objective(dat, fit$u_mat, fit$v_mat)
    res2 <- .evaluate_objective(dat, matrix(1, ncol = 1, nrow = 5),
                                matrix(1, ncol = 1, nrow = 5))
    res3 <- .evaluate_objective(dat, abs(matrix(rnorm(5), ncol = 1, nrow = 5)),
                                abs(matrix(rnorm(5), ncol = 1, nrow = 5)))

    res1 < res2 & res1 < res3
  })

  expect_true(all(bool_vec))

})



test_that("fit_factorization is appropriately minimized among rank 1 matrices", {
  trials <- 10
  n <- 20

  dat_list <- lapply(1:trials, function(x){
    set.seed(x)
    u_mat <- matrix(3*rnorm(n), ncol = 1)
    v_mat <- matrix(3*rnorm(n), ncol = 1)
    pred_mat <- u_mat %*% t(v_mat)

    dat <- pred_mat

    for(i in 1:n){
      for(j in 1:n){
        dat[i,j] <- abs(stats::rnorm(1, mean = pred_mat[i,j], sd = 0.5))
      }
    }

    attr(dat, "family") <- "gaussian"

    dat
  })

  fit_list <- lapply(1:trials, function(i){
    init <- initialization(dat_list[[i]], k = 1, family = "gaussian", max_val = 100)

    fit <- fit_factorization(dat_list[[i]], k = 1, u_mat = init$u_mat, v_mat = init$v_mat,
                             max_iter = 10, max_val = 100,
                             family = "gaussian")
  })

  # compare the grid of objective values to dataset
  error_mat <- sapply(1:trials, function(i){
    vec <- sapply(1:trials, function(j){
      .evaluate_objective(dat_list[[i]], fit_list[[j]]$u_mat, fit_list[[j]]$v_mat)
    })
    (vec - min(vec))/diff(range(vec))
  })

  expect_true(all(diag(error_mat) <= 1e-6))
})


test_that("fit_factorization is about the (roughly) same as SVD", {
  trials <- 10

  diff_val <- sapply(1:trials, function(x){
    set.seed(11*x)
    dat <- matrix(abs(rnorm(25, 2, 1)), nrow = 5, ncol = 5)
    attr(dat, "family") <- "gaussian"
    init <- initialization(dat, family = "gaussian", max_val = 100)

    fit <- fit_factorization(dat, k=2, u_mat = init$u_mat, v_mat = init$v_mat,
                             max_iter = 50, max_val = 100,
                             family = "gaussian")
    pred_mat1 <- fit$u_mat %*% t(fit$v_mat)

    svd_res <- svd(dat)
    pred_mat2 <- svd_res$u[,1:2]%*%diag(svd_res$d[1:2])%*%t(svd_res$v[,1:2])

    sum(abs(pred_mat1 - pred_mat2))
  })

  expect_true(sum(abs(diff_val)) <= 1e-6)
})





