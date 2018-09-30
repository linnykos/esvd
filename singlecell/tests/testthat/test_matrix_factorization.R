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

test_that(".gradient_vec satisfies the gradient definition", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
    u_vec <- abs(rnorm(2))
    u_vec2 <- abs(rnorm(2))
    v_mat <- abs(matrix(rnorm(8), nrow = 4, ncol = 2))

    i <- sample(1:10, 1)
    grad <- .gradient_vec(dat[i,], u_vec, v_mat)

    res <- .evaluate_objective_single(dat[i,], u_vec, v_mat)
    res2 <- .evaluate_objective_single(dat[i,], u_vec2, v_mat)

    res2 >= res + as.numeric(grad %*% (u_vec2 - u_vec)) - 1e-6
  })

  res2 - (res + as.numeric(grad %*% (u_vec2 - u_vec)))

  expect_true(all(bool_vec))
})


##################

## .evaluate_objective is correct

test_that(".evaluate_objective works", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
  v_mat <- abs(matrix(rnorm(8), nrow = 4, ncol = 2))

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
    v_mat <- abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    pred_mat <- u_mat %*% t(v_mat)
    dat_mat <- pred_mat

    for(i in 1:10){
      for(j in 1:4){
        dat_mat[i,j] <- stats::rnorm(1, pred_mat[i,j], pred_mat[i,j]/2)
      }
    }

    res <- .evaluate_objective(dat_mat, u_mat, v_mat)

    u_mat2 <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
    v_mat2 <- abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    res2 <- .evaluate_objective(dat_mat, u_mat2, v_mat2)

    c(res, res2)
  })

  expect_true(mean(avg_obj[1,]) < mean(avg_obj[2,]))
})

test_that(".evaluate_objective is equal to many .evaluate_objective_single", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
  v_mat <- abs(matrix(rnorm(8), nrow = 4, ncol = 2))

  res <- .evaluate_objective(dat, u_mat, v_mat)

  res2 <- sum(sapply(1:nrow(u_mat), function(x){
    .evaluate_objective_single(dat[x,], u_mat[x,], v_mat)
  }))

  expect_true(abs(res - res2) <= 1e-6)
})

test_that(".evaluate_objective gives sensible optimal", {
  set.seed(20)
  dat <- abs(matrix(rnorm(100, mean = 50, sd = 50/2), nrow = 10, ncol = 10))
  u_mat <- matrix(10, nrow = 10, ncol = 1)
  v_mat <- matrix(5, nrow = 10, ncol = 1)

  res <- .evaluate_objective(dat, u_mat, v_mat)

  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    u_mat2 <- abs(matrix(rnorm(10, mean = 10), nrow = 10, ncol = 1))
    v_mat2 <- abs(matrix(rnorm(10, mean = 5), nrow = 10, ncol = 1))
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
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA
  u_mat <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
  v_mat <- abs(matrix(rnorm(8), nrow = 4, ncol = 2))

  i <- 5
  res <- .evaluate_objective_single(dat[i,], u_mat[i,], v_mat)

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
    v_mat <- abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    pred_mat <- u_mat %*% t(v_mat)
    dat_mat <- pred_mat

    for(i in 1:10){
      for(j in 1:4){
        dat_mat[i,j] <- stats::rnorm(1, pred_mat[i,j], pred_mat[i,j]/2)
      }
    }

    i <- sample(1:10, 1)
    res <- .evaluate_objective(dat_mat[i,], u_mat[i,], v_mat)

    u_mat2 <- abs(matrix(rnorm(20), nrow = 10, ncol = 2))
    v_mat2 <- abs(matrix(rnorm(8), nrow = 4, ncol = 2))
    res2 <- .evaluate_objective(dat_mat[i,], u_mat2[i,], v_mat2)

    c(res, res2)
  })

  expect_true(mean(avg_obj[1,]) < mean(avg_obj[2,]))
})

################

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

test_that(".projection_l1 maintains greater than 0", {
  set.seed(15)
  current_vec <- rnorm(10)
  other_mat <- matrix(rnorm(50), ncol = 10)
  res <- .projection_l1(current_vec, other_mat)

  expect_true(all(other_mat %*% res >= -1e-6))
})

test_that(".projection_l1 can keep the current vector", {
  set.seed(20)
  current_vec <- abs(rnorm(10))
  other_mat <- abs(matrix(rnorm(50), ncol = 10))
  res <- .projection_l1(current_vec, other_mat)

  expect_true(sum(abs(res - current_vec)) < 1e-6)
})

test_that(".projection_l1 is actually a projection compared to the all 0 vector", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    current_vec <- abs(rnorm(10))
    other_mat <- abs(matrix(rnorm(50), ncol = 10))
    new_vec <- .projection_l1(current_vec, other_mat)

    sum(abs(current_vec - new_vec)) <= sum(abs(current_vec))
  })

  expect_true(all(bool_vec))
})

#############

## .initialization is correct

test_that(".initialization works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA

  res <- .initialization(dat)

  expect_true(is.list(res))
  expect_true(ncol(res$u_mat) == ncol(res$v_mat))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
})

test_that(".initialization actually gives positive predictions", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA

  res <- .initialization(dat)
  pred_mat <- res$u_mat %*% t(res$v_mat)

  expect_true(all(pred_mat[which(!is.na(dat))] > -1e6))
})

#########################

## .backtrack_linesearch is correct

test_that(".backtrack_linesearch works", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA

  res <- .initialization(dat)
  u_mat <- res$u_mat
  v_mat <- res$v_mat
  i <- 1
  grad_vec <- .gradient_vec(dat[i,], u_mat[i,], v_mat)

  res <- .backtrack_linesearch(dat[i,], u_mat[i,], v_mat, grad_vec)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
})

test_that(".backtrack_linesearch actually keeps the positive constraint", {
  trials <- 50

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)

    dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
    dat[sample(1:prod(dim(dat)), 10)] <- NA

    res <- .initialization(dat)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    i <- sample(1:10, 1)

    if(any(!is.na(dat[i,]))){
      grad_vec <- .gradient_vec(dat[i,], u_mat[i,], v_mat)

      res <- .backtrack_linesearch(dat[i,], u_mat[i,], v_mat, grad_vec)

      u_new <- u_mat[i,] - res*grad_vec
      all(v_mat %*% u_new >= -1e-6)
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
    dat[sample(1:prod(dim(dat)), 10)] <- NA

    res <- .initialization(dat)
    u_mat <- res$u_mat
    v_mat <- res$v_mat

    i <- sample(1:10, 1)
    if(any(!is.na(dat[i,]))){
      obj1 <- .evaluate_objective_single(dat[i,], u_mat[i,], v_mat)
      grad_vec <- .gradient_vec(dat[i,], u_mat[i,], v_mat)

      res <- .backtrack_linesearch(dat[i,], u_mat[i,], v_mat, grad_vec)

      u_new <- u_mat[i,] - res*grad_vec
      obj2 <- .evaluate_objective_single(dat[i,], u_new, v_mat)

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
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA

  res <- .initialization(dat)
  u_mat <- res$u_mat
  v_mat <- res$v_mat
  i <- 1

  res <- .optimize_row(dat[i,], u_mat[i,], v_mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == length(u_mat[i,]))
})

test_that(".optimize_row actually lowers the objective", {
  trials <- 25

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
    dat[sample(1:prod(dim(dat)), 10)] <- NA

    res <- .initialization(dat)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    i <- sample(1:10, 1)

    if(any(!is.na(dat[i,]))){
      u_new <- .optimize_row(dat[i,], u_mat[i,], v_mat)
      obj1 <- .evaluate_objective_single(dat[i,], u_mat[i,], v_mat)
      obj2 <- .evaluate_objective_single(dat[i,], u_new, v_mat)

      obj2 <= obj1 + 1e-6
    } else {TRUE}
  })

  expect_true(all(bool_vec))
})

test_that(".optimize_row works the other way", {
  trials <- 25

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
    dat[sample(1:prod(dim(dat)), 10)] <- NA

    res <- .initialization(dat)
    u_mat <- res$u_mat
    v_mat <- res$v_mat
    j <- sample(1:4, 1)

    if(any(!is.na(dat[,j]))){
      v_new <- .optimize_row(dat[,j], v_mat[j,], u_mat)
      obj1 <- .evaluate_objective_single(dat[,j], v_mat[j,], u_mat)
      obj2 <- .evaluate_objective_single(dat[,j], v_new, u_mat)

      obj2 <= obj1 + 1e-6
    } else {TRUE}
  })

  expect_true(all(bool_vec))
})

##################

## .optimize_mat is correct

test_that(".optimize_mat works", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  dat[sample(1:prod(dim(dat)), 10)] <- NA

  res <- .initialization(dat)
  u_mat <- res$u_mat
  v_mat <- res$v_mat

  res <- .optimize_mat(dat, u_mat, v_mat)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(u_mat)))
})

test_that(".optimize_mat keeps the positive constraint", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
    dat[sample(1:prod(dim(dat)), 10)] <- NA
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

    all(pred_mat[idx] >= -1e-6)
  })

  expect_true(all(bool_vec))
})


test_that(".optimize_mat lowers the objective value", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
    dat[sample(1:prod(dim(dat)), 10)] <- NA
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

## .fit_gaussian_factorization is correct

test_that(".fit_gaussian_factorization works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(20), nrow = 5, ncol = 4))

  res <- .fit_gaussian_factorization(dat, k = 2)

  expect_true(is.list(res))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
  expect_true(ncol(res$u_mat) == 2)
  expect_true(ncol(res$v_mat) == 2)
})

test_that(".fit_gaussian_factorization works with missing values", {
  set.seed(5)
  dat <- abs(matrix(rnorm(20), nrow = 5, ncol = 4))
  dat[sample(1:prod(dim(dat)), 5)] <- NA

  res <- .fit_gaussian_factorization(dat)

  expect_true(is.list(res))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
  expect_true(ncol(res$u_mat) == 2)
  expect_true(ncol(res$v_mat) == 2)
})

test_that(".fit_gaussian_factorization preserves the positive entries", {
  trials <- 10

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x*10)
    dat <- abs(matrix(rnorm(20), nrow = 5, ncol = 4))
    dat[sample(1:prod(dim(dat)), 5)] <- NA

    res <- .fit_gaussian_factorization(dat)

    pred_mat <- res$u_mat %*% t(res$v_mat)

    all(pred_mat[which(!is.na(dat))] >= -1e-6)
  })

  expect_true(all(bool_vec))

})
