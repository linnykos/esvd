context("Test slingshot 3D plots")

## .remove_duplicate_rows is correct

test_that(".remove_duplicate_rows works", {
  set.seed(10)
  mat <- matrix(1:100, nrow = 20)
  duplicate_vec <- sapply(1:20, function(x){sample(1:5,1)})
  mat2 <- do.call(rbind, lapply(1:20, function(j){
    matrix(rep(mat[j,], times = duplicate_vec[j]), nrow = duplicate_vec[j], byrow = T)
  }))

  res <- .remove_duplicate_rows(mat2)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(mat)))
  expect_true(sum(abs(res - mat)) <= 1e-5)
})

##############################

## .find_adjacent_directions is correct

test_that(".find_adjacent_directions works", {
  set.seed(10)
  dat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  idx <- 5
  res <- .find_adjacent_directions(dat, idx)

  expect_true(length(res) == ncol(dat))
  expect_true(is.numeric(res))
})

test_that(".find_adjacent_directions can work for idx = 1", {
  set.seed(10)
  dat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  idx <- 1
  res <- .find_adjacent_directions(dat, idx)

  expect_true(length(res) == ncol(dat))
  expect_true(is.numeric(res))
})

test_that(".find_adjacent_directions can work for idx = nrow(dat)", {
  set.seed(10)
  dat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  idx <- nrow(dat)
  res <- .find_adjacent_directions(dat, idx)

  expect_true(length(res) == ncol(dat))
  expect_true(is.numeric(res))
})

test_that(".find_adjacent_directions produces unit vectors", {
  trials <- 20

  bool_vec <- sapply(1:trials, function(trial){
    set.seed(trial)
    dat <- matrix(rnorm(100), nrow = 20, ncol = 5)
    idx <- nrow(dat)
    res <- .find_adjacent_directions(dat, idx)

    abs(.l2norm(res) - 1) <= 1e-5
  })

  expect_true(all(bool_vec))
})

test_that(".find_adjacent_directions gives the same direction for a set of points along a line", {
  mat <- matrix(1:100, nrow = 20)
  res1 <- .find_adjacent_directions(mat, idx = 5)
  res2 <- .find_adjacent_directions(mat, idx = 15)

  expect_true(sum(abs(res1 - res2)) <= 1e-5)
})

##################################

## .projection_matrix is correct

test_that(".projection_matrix works", {
  set.seed(10)
  vec <- rnorm(20)
  mat <- matrix(rnorm(100), nrow = 20, ncol = 5)

  res <- .projection_matrix(vec, mat)
  expect_true(is.numeric(res))
  expect_true(length(res) == length(vec))
})

test_that(".projection_matrix doesn't change when applied twice", {
  set.seed(5)
  vec <- rnorm(20)
  mat <- matrix(rnorm(100), nrow = 20, ncol = 5)

  res <- .projection_matrix(vec, mat)
  res2 <- .projection_matrix(res, mat)

  expect_true(sum(abs(res - res2)) <= 1e-5)
})

##################################

## .find_basis_vectors is correct

test_that(".find_basis_vectors works", {
  res <- .find_basis_vectors(c(1:3))

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == c("vec1", "vec2")))
})

test_that(".find_basis_vectors finds vectors that are perpendicular", {
  vec <- c(1:3)
  res <- .find_basis_vectors(vec)

  val1 <- abs(vec%*%res$vec1)
  val2 <- abs(vec%*%res$vec2)
  val3 <- abs(res$vec1%*%res$vec2)

  expect_true(val1 + val2 + val3 <= 1e-5)
})

test_that(".find_basis_vectors finds unit vectors", {
  vec <- c(1:3)
  res <- .find_basis_vectors(vec)

  val1 <- .l2norm(res$vec1)
  val2 <- .l2norm(res$vec2)

  expect_true(abs(val1 - 1) <= 1e-5)
  expect_true(abs(val2 - 1) <= 1e-5)
})

########################################

## .construct_3d_circle is correct

test_that(".construct_3d_circle works",{
  dat_vec <- c(1,2,3)
  radius <- 1
  basis_res <- .find_basis_vectors(c(1:3))
  res <- .construct_3d_circle(dat_vec, radius, basis_vec1 = basis_res$vec1,
                              basis_vec2 = basis_res$vec2, len = 20)

  expect_true(all(dim(res) == c(20, 3)))
  expect_true(is.matrix(res))
})

test_that(".construct_3d_circle makes a circle",{
  dat_vec <- c(1,2,3)
  radius <- 1
  basis_res <- .find_basis_vectors(c(1:3))
  res <- .construct_3d_circle(dat_vec, radius, basis_vec1 = basis_res$vec1,
                              basis_vec2 = basis_res$vec2, len = 20)

  bool_vec <- sapply(1:nrow(res), function(i){
    abs(.l2norm(dat_vec - res[i,]) - radius) <= 1e-5
  })

  expect_true(all(bool_vec))
})

test_that(".construct_3d_circle is in the plane of the basis vectors",{
  dat_vec <- c(1,2,3)
  radius <- 1
  basis_res <- .find_basis_vectors(c(1:3))
  res <- .construct_3d_circle(dat_vec, radius, basis_vec1 = basis_res$vec1,
                              basis_vec2 = basis_res$vec2, len = 20)

  basis_mat <- cbind(basis_res$vec1, basis_res$vec2, c(1:3)/.l2norm(c(1:3)))

  bool_vec <- sapply(1:nrow(res), function(i){
    tmp <- solve(basis_mat, res[i,] - dat_vec)
    abs(tmp[3]) <= 1e-5
  })

  expect_true(all(bool_vec))
})

#####################################3

## .construct_all_circles is correct

test_that(".construct_all_circles works", {
  mat <- matrix(1:90, nrow = 30, ncol = 3)
  radius <- 1
  res <- .construct_all_circles(mat, radius)

  expect_true(length(res) == nrow(mat))
  expect_true(is.list(res))
  expect_true(all(sapply(res, nrow) == 20))
})

##########################################3

## .find_correct_orientation is correct

test_that(".find_correct_orientation works", {
  mat <- matrix(1:90, nrow = 30, ncol = 3)
  radius <- 1
  circle_list <- .construct_all_circles(mat, radius)
  mat1 <- circle_list[[1]]
  mat2 <- circle_list[[2]]

  res <- .find_correct_orientation(mat1, mat2)

  expect_true(all(dim(res) == dim(mat2)))
  expect_true(sum(res) == sum(mat2))
})

############################################

## .subset_indices is correct

test_that(".subset_indices works", {
  res <- .subset_indices(10, 2, 6, direction_forward = T)

  expect_true(is.list(res))
  expect_true(all(sort(res$to) == 1:10))
  expect_true(all(sort(res$from) == 1:10))
})

#############################################

## construct_3d_tube is correct

test_that("construct_3d_tube works", {
  mat <- matrix(1:90, nrow = 30, ncol = 3)
  radius <- 1
  res <- construct_3d_tube(mat, radius, len = 20)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == c("x_mat", "y_mat", "z_mat")))
  expect_true(all(dim(res$x_mat) == c(20, nrow(mat))))
  expect_true(all(dim(res$y_mat) == c(20, nrow(mat))))
  expect_true(all(dim(res$z_mat) == c(20, nrow(mat))))
})
