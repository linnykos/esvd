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
  expect_true(all(names(res$pcurve_list) == paste0("Lineage", 1:length(res$pcurve_list))))
})

########################

## .smoother_func is correct

test_that(".smoother_func works", {
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
  pcurve_list <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list

  sample_idx <- .determine_idx_lineage(lineages[[1]], cluster_mat)
  res <- .smoother_func(pcurve_list[[1]]$lambda, dat[sample_idx,], b = 1)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(dat[sample_idx,])))
})

########################

## .construct_average_curve is correct

test_that(".construct_average_curve works", {
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
  pcurve_list <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list

  res <- .construct_average_curve(pcurve_list, dat)

  expect_true(class(res) == "principal_curve")
})

###############################

## .refine_curve_fit is correct

test_that(".refine_curve_fit works", {
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
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("pcurve_list", "D")))
  expect_true(all(dim(res$D) == c(100, length(lineages))))
  expect_true(is.list(res$pcurve_list))
  expect_true(all(sapply(res$pcurve_list, class) == "principal_curve"))
})

##############################

## .percent_shrinkage is correct

test_that(".percent_shrinkage works", {
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
  pcurve <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list[[1]]
  common_idx <- which(W[,1] == 1)

  res <- .percent_shrinkage(pcurve, common_idx)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == length(pcurve$lambda))
})


test_that(".percent_shrinkage computes a shrinkage for all indices, even outliers ", {
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
  pcurve <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list[[1]]
  pcurve$lambda[1:5] <- 10*max(pcurve$lambda)
  common_idx <- which(W[,1] == 1)

  res <- .percent_shrinkage(pcurve, common_idx)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == length(pcurve$lambda))
  expect_true(all(res[1:5] == 0))
})

test_that(".percent_shrinkage respects the order of lambda", {
  set.seed(20)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, knn = NA, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
  pcurve <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list[[1]]
  common_idx <- which(W[,1] == 1)

  res <- .percent_shrinkage(pcurve, common_idx)

  res <- res[order(pcurve$lambda, decreasing = F)]
  vec <- diff(res)
  expect_true(all(vec <= 0))
})

###########################

## .shrink_to_avg is correct

test_that(".shrink_to_avg works", {
  set.seed(20)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, knn = NA, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)


  pcurve_list <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list
  avg_curve <- .construct_average_curve(pcurve_list, dat)

  pcurve <- pcurve_list[[1]]
  common_idx <- which(W[,1] == 1)
  pct_shrink <- .percent_shrinkage(pcurve, common_idx)

  res <- .shrink_to_avg(pcurve, avg_curve, pct_shrink, dat)

  expect_true(class(res) == "principal_curve")
  expect_true(all(res$W == pcurve$W))

})
