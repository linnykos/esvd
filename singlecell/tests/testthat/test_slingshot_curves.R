context("Test slingshot curves")

## .initialize_weight_matrix is correct

test_that(".initialize_weight_matrix works", {
  set.seed(10)
  cluster_labels <- sample(1:10, 200, replace = T)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, knn = NA, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)

  res <- .initialize_weight_matrix(cluster_mat, lineages)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(200, length(lineages))))
  expect_true(all(sort(unique(as.numeric(res)))== c(0,1)))
})

###############

## .initialize_curve_hierarchy is correct

test_that(".initialize_curve_hierarchy works", {
  lineages <- list(Lineage1 = c(1,2,3), Lineage2 = c(1,2,4))
  cluster_vec <- 1:4

  res <- .initialize_curve_hierarchy(lineages, cluster_vec)

  expect_true(is.list(res))
})

test_that(".initialize_curve_hierarchy works for a more complicated case", {
  set.seed(10)
  cluster_labels <- sample(1:10, 200, replace = T)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, knn = NA, starting_cluster = 1)
  cluster_vec <- 1:10

  res <- .initialize_curve_hierarchy(lineages, cluster_vec)

  expect_true(is.list(res))
})

#################

## .initial_curve_fit is correct

test_that(".initial_curve_fit works", {
  set.seed(10)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, knn = NA, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)

  res <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)

  expect_true(is.list(res))
  expect_true(all(names(res) == c("pcurve_list", "D")))
  expect_true(all(sapply(res$pcurve_list, class) == "principal_curve"))
  expect_true(all(dim(res$D) == c(100, length(lineages))))
})
