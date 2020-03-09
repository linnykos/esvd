context("Test missing values")

test_that("construct_missing_values works", {
  set.seed(10)
  res <- construct_missing_values(10,50)
  expect_true(all(res %% 1 == 0))
  expect_true(all(res >= 1))
  expect_true(all(res <= 10*50))
})

test_that("construct_missing_values is indexing correctly", {
  set.seed(10)
  n <- 100; d<- 500
  res <- construct_missing_values(n, d)

  set.seed(10)
  missing_mat <- rbind(do.call(rbind, (lapply(1:n, function(x){
    cbind(x, sample(1:d, 4))
  }))), do.call(rbind, (lapply(1:d, function(x){
    cbind(sample(1:n, 4), x)
  }))))

  mat <- matrix(1, nrow = n, ncol = d)
  for(i in 1:nrow(missing_mat)){
    mat[missing_mat[i,1], missing_mat[i,2]] <- NA
  }
  res2 <- which(is.na(mat))

  expect_true(all(sort(res) == sort(res2)))
})

test_that("construct_missing_values respects num_val", {
  n <- 100; d<- 500
  set.seed(20)
  res <- construct_missing_values(n, d, 1)

  set.seed(20)
  res2 <- construct_missing_values(n, d, 2)

  expect_true(length(res) < length(res2))
})
