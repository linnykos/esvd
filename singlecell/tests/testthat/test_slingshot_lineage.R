context("Test slingshot lineages")

## .construct_cluster_matrix is correct

test_that(".construct_cluster_matrix works", {
  set.seed(10)
  res <- .construct_cluster_matrix(sample(1:10, 200, replace = T))

  expect_true(all(dim(res) == c(200,10)))
  expect_true(all(sort(unique(as.numeric(res))) == c(0,1)))
})

###########

## .compute_cluster_center is correct

test_that(".compute_cluster_center works", {
  set.seed(10)
  cluster_mat <- .construct_cluster_matrix(sample(1:10, 200, replace = T))
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))

  res <- .compute_cluster_center(dat, cluster_mat)

  expect_true(all(dim(res) == c(10,5)))
})

###############

## .construct_knn_graph is correct

test_that(".construct_knn_graph works", {
  set.seed(10)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))

  res <- .construct_knn_graph(dat, 5)

  expect_true(class(res) == "igraph")
})

test_that(".construct_knn_graph is actually constructing the graph based on min distance", {
  dat <- t(sapply(1:10, function(x){
    rep(x, 5)
  }))

  res <- .construct_knn_graph(dat, 2)
  adj <- as.matrix(igraph::as_adjacency_matrix(res))

  expect_true(all(diag(adj) == 0))
  expect_true(adj[1,2] == 1)
})
