context("Test dropout")

## dropout_create is correct

test_that("dropout_create works", {
  fun <- dropout_create(1,1)
  vec <- sapply(seq(-10,10,length.out=100), fun)

  expect_true(length(vec) == 100)
  expect_true(is.numeric(vec))
  expect_true(all(diff(vec) >= 0))
})

############

## .form_adj_nonzero is correct

test_that(".form_adj_nonzero works", {
  set.seed(10)
  dat <- matrix(1:200,10,20)
  dat[sample(1:200, 100)] <- 0

  res <- .form_adj_nonzero(dat)

  n <- nrow(dat)
  d <- ncol(dat)
  expect_true(all(dim(res) == rep(n+d, 2)))
  expect_true(all(sum(dat) == sum(dat*res[1:n,(n+1:d)])))
  expect_true(all(sort(unique(as.numeric(res))) == c(0,1)))
})

#################

## .largest_average_degree_subgraph is correct

test_that(".largest_average_degree_subgraph works", {
  set.seed(10)
  adj <- matrix(0, 30, 30)
  adj[sample(1:30^2, 400)] <- 1
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  igraph::V(g)$name <- as.character(1:30)

  res <- .largest_average_degree_subgraph(g, 1)

  expect_true(length(res) <= 30)
  expect_true(all(is.character(res)))
})

test_that(".largest_average_degree_subgraph works even if graph is not labeled", {
  set.seed(10)
  adj <- matrix(0, 30, 30)
  adj[sample(1:30^2, 400)] <- 1
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")

  res <- .largest_average_degree_subgraph(g, 1)

  expect_true(length(res) <= 30)
  expect_true(all(is.character(res)))
})

###################

## .funk_svd_prediction is correct

test_that(".funk_svd_prediction works", {
  set.seed(10)
  dat <- matrix(1:200,10,20)
  dat[sample(1:200, 100)] <- 0

  res <- .funk_svd_prediction(dat, lambda = 0.0001)

  expect_true(all(dim(res) == dim(dat)))
  expect_true(is.matrix(res))
  expect_true(!any(is.na(res)))
  expect_true(!any(is.nan(res)))
})

###################

## .logistic_regression is correct

test_that(".logistic_regression works", {
  set.seed(10)
  observed_dat <- rnorm(100)
  observed_dat[observed_dat < 0] <- 0
  pred_dat <- rnorm(100)

  res <- .logistic_regression(observed_dat, pred_dat, threshold_quant = 0.1)

  expect_true(length(res) == 2)
  expect_true(is.numeric(res))
  expect_true(all(sort(names(res)) == c("Intercept", "Slope")))
})

####################

## estimate_dropout is correct

test_that("estimate_dropout works", {
  set.seed(10)
  dat <- matrix(1:200,10,20)
  dat[sample(1:200, 100)] <- 0

  res <- estimate_dropout(dat, lambda = 0.0001)

  expect_true(length(res) == 2)
  expect_true(is.numeric(res))
  expect_true(!any(is.na(res)))
  expect_true(!any(is.nan(res)))
})

test_that("estimate_drop will crash if lambda is too large", {
  set.seed(10)
  dat <- matrix(1:200,10,20)
  dat[sample(1:200, 100)] <- 0

  expect_error(estimate_dropout(dat, lambda = 0.01))
})
