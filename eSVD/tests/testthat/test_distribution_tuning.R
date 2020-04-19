context("Test tuning the distribution")

## plot_prediction_against_observed is correct

test_that("plot_prediction_against_observed works", {
  set.seed(10)
  dat <- matrix(rnbinom(40, size = 10, prob = 0.5), nrow = 5, ncol = 5)
  missing_idx <- construct_missing_values(n = nrow(dat), p = ncol(dat), num_val = 1)
  dat_NA <- dat
  dat_NA[missing_idx] <- NA

  init <- initialization(dat_NA, family = "neg_binom", max_val = 100, scalar = 10)
  fit <- fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                           max_iter = 10, max_val = 100,
                           family = "neg_binom", scalar = 10)

  res <- plot_prediction_against_observed(dat, nat_mat_list = list(fit$u_mat %*% t(fit$v_mat)),
                                          family = "neg_binom", missing_idx_list = list(missing_idx),
                                          scalar = 10, plot = F)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(sort(names(res)) == sort(c("angle_val", "bool"))))
  expect_true(is.numeric(res$angle_val))
  expect_true(length(res$angle_val) == 1)
  expect_true(res$angle_val >= 0)
  expect_true(res$angle_val <= 180)
  expect_true(length(res$bool) == 1)
  expect_true(is.logical(res$bool))
})

##############

## tuning_select_scalar is correct

test_that("tuning_select_scalar works", {
  set.seed(10)
  dat <- matrix(rnbinom(200, size = 10, prob = 0.5), nrow = 20, ncol = 10)
  missing_idx <- construct_missing_values(n = nrow(dat), p = ncol(dat), num_val = 1)
  dat_NA <- dat
  dat_NA[missing_idx] <- NA

  scalar_vec <- c(2, 10, 100)

  fit_list <- lapply(scalar_vec, function(scalar){
    init <- initialization(dat_NA, family = "neg_binom", max_val = 100, scalar = scalar, k = 1)
    fit <- fit_factorization(dat_NA, u_mat = init$u_mat, v_mat = init$v_mat,
                             max_iter = 10, max_val = 100, k = 1,
                             family = "neg_binom", scalar = scalar)
  })

  nat_mat_list_list <- lapply(1:length(scalar_vec), function(i){
    list(fit_list[[i]]$u_mat %*% t(fit_list[[i]]$v_mat))
  })

  res <- tuning_select_scalar(dat, nat_mat_list = nat_mat_list_list,
                              family = "neg_binom", missing_idx_list = list(missing_idx),
                              scalar_vec = scalar_vec)


})

#################################

## .compute_principal_angle is correct

test_that(".compute_principal_angle works", {
  set.seed(10)
  tmp_mat <- cbind(1:10, 1:10+rnorm(10))
  res <- .compute_principal_angle(tmp_mat)

  expect_true(res >= 0)
  expect_true(res <= 180)
})

test_that(".compute_principal_angle can handle when the angle is 45", {
  set.seed(10)
  tmp_mat <- cbind(1:10, 1:10)
  res <- .compute_principal_angle(tmp_mat)

  expect_true(abs(res - 45) <= 1e-6)
})

test_that(".compute_principal_angle can handle when the angle is 45, even when data is reflected", {
  set.seed(10)
  tmp_mat <- cbind(-(1:10), -(1:10))
  res <- .compute_principal_angle(tmp_mat)

  expect_true(abs(res - 45) <= 1e-6)
})

test_that(".compute_principal_angle can handle when the angle is slightly above 45", {
  set.seed(10)
  tmp_mat <- cbind((1:10)+0.1*(1:10), 1:10)
  res <- .compute_principal_angle(tmp_mat)

  expect_true(res > 45)
})


test_that(".compute_principal_angle can handle when the angle is 0", {
  set.seed(10)
  tmp_mat <- cbind(rep(0, 10), 0:9)
  res <- .compute_principal_angle(tmp_mat)

  expect_true(abs(res - 0) <= 1e-6 || abs(res - 180) <= 1e-6)
})

test_that(".compute_principal_angle can handle when the angle is 90", {
  set.seed(10)
  tmp_mat <- cbind(0:9, rep(0, 10))
  res <- .compute_principal_angle(tmp_mat)

  expect_true(abs(res - 90) <= 1e-6)
})

test_that(".compute_principal_angle can handle negative angles", {
  set.seed(10)
  tmp_mat <- cbind(-(0:10), 0:10)
  res <- .compute_principal_angle(tmp_mat)
  expect_true(abs(res - 135) <= 1e-6)

  set.seed(10)
  tmp_mat <- cbind((0:10), -c(0:10))
  res <- .compute_principal_angle(tmp_mat)

  expect_true(abs(res - 135) <= 1e-6)
})

#####################3

## .within_prediction_region is correct

test_that(".within_prediction_region works", {
  res <- .within_prediction_region(50, family = "neg_binom", width = 0.8, scalar = 10, angle_val = 46)

  expect_true(is.list(res))
  expect_true(length(res) == 4)
  expect_true(all(sort(names(res)) == sort(c("seq_vec", "interval_mat", "principal_line", "bool"))))
  expect_true(length(res$seq_vec) == length(res$principal_line))
  expect_true(length(res$seq_vec) == ncol(res$interval_mat))
  expect_true(is.logical(res$bool))
})

test_that(".within_prediction_region outputs the correct answer in an easy 45 degree setting", {
  res <- .within_prediction_region(50, family = "neg_binom", width = 0.8, scalar = 10, angle_val = 45)

  expect_true(res$bool)
})

test_that(".within_prediction_region will still work with negative angles", {
  res <- .within_prediction_region(50, family = "neg_binom", width = 0.8, scalar = 10, angle_val = 135)

  expect_true(all(res$principal_line <= 0))
  expect_true(!res$bool)
})


test_that(".within_prediction_region outputs the correct answer in an easy 47 degree setting", {
  res <- .within_prediction_region(50, family = "neg_binom", width = 0.8, scalar = 10, angle_val = 45)

  expect_true(res$bool)
})

test_that(".within_prediction_region correctly works within the workflow", {
  set.seed(10)
  dat <- matrix(rnbinom(40, size = 10, prob = 0.5), nrow = 5, ncol = 5)
  missing_idx <- construct_missing_values(n = nrow(dat), p = ncol(dat), num_val = 1)
  dat_NA <- dat
  dat_NA[missing_idx] <- NA

  init <- initialization(dat_NA, family = "neg_binom", max_val = 100, scalar = 10)
  fit <- fit_factorization(dat, u_mat = init$u_mat, v_mat = init$v_mat,
                           max_iter = 10, max_val = 100,
                           family = "neg_binom", scalar = 10)

  nat_mat_list = list(fit$u_mat %*% t(fit$v_mat))
  family = "neg_binom"
  missing_idx_list = list(missing_idx)
  scalar = 10
  width = 0.8
  plot = F
  tol = 0.95
  max_points = 500000

  stopifnot(length(nat_mat_list) == length(missing_idx_list))

  pred_mat_list <- lapply(nat_mat_list, function(nat_mat){
    compute_mean(nat_mat, family = family, scalar = scalar)
  })

  tmp_list <- lapply(1:length(nat_mat_list), function(i){
    cbind(dat[missing_idx_list[[i]]], pred_mat_list[[i]][missing_idx_list[[i]]])
  })

  angle_vec <- sapply(tmp_list, .compute_principal_angle)
  angle_val <- mean(angle_vec)

  tmp_mat <- do.call(rbind, tmp_list)
  colnames(tmp_mat) <- c("observed_val", "predicted_val")
  res <- .within_prediction_region(max(tmp_mat[,"predicted_val"]), family = family, width = width,
                                   scalar = scalar, angle_val = angle_val, tol = tol)

  expect_true(res$bool)
})
