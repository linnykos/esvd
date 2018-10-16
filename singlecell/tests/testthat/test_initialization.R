context("Test initialization")

## .find_neighbors_impute is correct

test_that(".find_neighbors_impute works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))

  res <- .find_neighbors_impute(dat, Kcluster = 2)

  expect_true(is.vector(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == nrow(dat))
  expect_true(length(unique(res)) == 2)
})

test_that(".find_neighbors_impute can find the right neighbors", {
  set.seed(10)
  dat <- rbind(MASS::mvrnorm(20, rep(0, 10), diag(10)),
               MASS::mvrnorm(20, rep(2, 10), 2*diag(10)))

  res <- .find_neighbors_impute(dat, Kcluster = 2)

  expect_true(length(unique(res[1:20])) == 1)
  expect_true(length(unique(res[21:40])) == 1)
})

############################

## .nnls_impute is correct

test_that(".nnls_impute works", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 4, ncol = 10))

  res <- .nnls_impute(dat[1,], dat[-1,], sample(1:10, 5))

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == ncol(dat))
})

test_that(".nnls_impute works for larger matrices", {
  set.seed(20)
  dat <- abs(matrix(rnorm(120), nrow = 4, ncol = 30))

  res <- .nnls_impute(dat[1,], dat[-1,], sample(1:30, 20))

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == ncol(dat))
})

test_that(".nnls_impute preserves the values that are kept", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 4, ncol = 10))
  vec <- rep(NA, 10)
  beta_vec <- rep(0.1, 4)
  vec <- as.numeric(beta_vec%*%dat)

  idx <- sample(1:10, 5)
  res <- .nnls_impute(vec, dat, idx)

  expect_true(all(vec[idx] == res[idx]))
})

test_that(".nnls_impute can reduce to linear regression", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 4, ncol = 10))
  vec <- rep(NA, 10)
  beta_vec <- rep(0.1, 4)
  vec <- as.numeric(beta_vec%*%dat)

  idx <- sample(1:10, 5)
  res <- .nnls_impute(vec, dat, idx)

  tmp <- as.data.frame(cbind(t(dat[,idx]), vec[idx]))
  colnames(tmp) <- c(paste0("X", 1:nrow(dat)), "Y")

  beta <- lm(Y~.-1, data = tmp)
  pred_val <- stats::coef(beta) %*% dat[,-idx]

  expect_true(sum(abs(sort(res[-idx]) == sort(pred_val))) <= 1e-6)
})

test_that(".nnls_impute can work when only two known genes are included", {
  set.seed(20)
  dat <- abs(matrix(rnorm(40), nrow = 4, ncol = 10))
  vec <- abs(rnorm(10))

  res <- .nnls_impute(vec, dat, 1:2, max_time = 5)

  expect_true(is.vector(res))
})

test_that(".nnls_impute does not get stuck in a loop", {
  set.seed(10)
  dat <- abs(rbind(MASS::mvrnorm(20, rep(0, 10), diag(10)),
                   MASS::mvrnorm(20, rep(2, 10), 2*diag(10))))
  drop_idx <- sample(1:prod(dim(dat)), 10)

  Kcluster <- 2
  neigh_vec <- .find_neighbors_impute(dat, Kcluster = Kcluster)
  neigh_list <- lapply(1:Kcluster, function(k){which(neigh_vec == k)})

  dat2 <- dat
  dat2[drop_idx] <- NA

  k <- neigh_list[[1]]
  i <- 23

  keep_idx <- which(!is.na(dat2[i,]))
  res <- .nnls_impute(dat[i,], dat[setdiff(k, i),,drop = F], keep_idx,
                      max_time = 5)

  expect_true(is.numeric(res))
})

######################

## .scImpute is correct

test_that(".scImpute works", {
  set.seed(10)
  dat <- abs(rbind(MASS::mvrnorm(5, rep(0, 5), diag(5)),
               MASS::mvrnorm(5, rep(10, 5), 2*diag(5))))
  dat2 <- dat
  for(i in 1:nrow(dat2)){
    dat2[i, sample(1:5, 1)] <- NA
  }
  drop_idx <- which(is.na(dat2))

  res <- .scImpute(dat, drop_idx = drop_idx, Kcluster = 2, min_size = 3,
                   max_time = 5)

  expect_true(nrow(res) == nrow(dat))
  expect_true(ncol(res) == ncol(res))
  expect_true(sum(abs(res[-drop_idx] - dat[-drop_idx])) <= 1e-6)
})



########################

## .initialization is correct

test_that(".initialization works", {
  set.seed(10)
  dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))

  res <- .initialization(dat)

  expect_true(is.list(res))
  expect_true(ncol(res$u_mat) == ncol(res$v_mat))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
})

test_that(".initialization works off of an imputed matrix", {
  set.seed(10)
  dat <- abs(rbind(MASS::mvrnorm(5, rep(0, 5), diag(5)),
                   MASS::mvrnorm(5, rep(10, 5), 2*diag(5))))
  dat2 <- dat
  for(i in 1:nrow(dat2)){
    dat2[i, sample(1:5, 1)] <- NA
  }
  drop_idx <- which(is.na(dat2))
  dat2 <- .scImpute(dat, drop_idx, Kcluster = 2, min_size = 3, max_time = 5)

  res <- .initialization(dat2)

  expect_true(is.list(res))
  expect_true(ncol(res$u_mat) == ncol(res$v_mat))
  expect_true(nrow(res$u_mat) == nrow(dat))
  expect_true(nrow(res$v_mat) == ncol(dat))
})

test_that(".initialization gives negative predictions", {
  set.seed(10)
  dat <- abs(rbind(MASS::mvrnorm(5, rep(0, 5), diag(5)),
                   MASS::mvrnorm(5, rep(10, 5), 2*diag(5))))
  dat2 <- dat
  for(i in 1:nrow(dat2)){
    dat2[i, sample(1:5, 1)] <- NA
  }
  drop_idx <- which(is.na(dat2))
  dat2 <- .scImpute(dat, drop_idx, Kcluster = 2, min_size = 3, max_time = 5)

  res <- .initialization(dat2)
  pred_mat <- res$u_mat %*% t(res$v_mat)

  expect_true(all(pred_mat[which(!is.na(dat))] <= 1e-6))
})

