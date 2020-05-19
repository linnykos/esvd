context("Test purity")

## .construct_neighborhood_graph is correct

test_that(".construct_neighborhood_graph works", {
  mat <- cbind(1:10, 1:10)
  neighborhood_size <- 3
  res <- .construct_neighborhood_graph(mat, neighborhood_size)

  expect_true(class(res) == "igraph")
})

test_that(".construct_neighborhood_graph removes duplicate edges", {
  set.seed(10)
  mat <- cbind(1:10, 1:10)
  order_vec <- sample(1:10)
  mat <- mat[order_vec,]
  neighborhood_size <- 3
  res <- .construct_neighborhood_graph(mat, neighborhood_size)

  adj_mat <- as.matrix(igraph::as_adjacency_matrix(res))

  expect_true(all(adj_mat %in% c(0,1)))
})

test_that(".construct_neighborhood_graph forms the correct graph", {
  set.seed(10)
  mat <- cbind(1:10, 1:10)
  order_vec <- sample(1:10)
  mat <- mat[order_vec,]
  neighborhood_size <- 2
  res <- .construct_neighborhood_graph(mat, neighborhood_size)

  edge_mat <- igraph::as_edgelist(res)
  edge_mat <- apply(edge_mat, 2, function(x){order_vec[x]})

  bool_vec <- apply(edge_mat, 1, function(x){
    abs(diff(x)) <= neighborhood_size
  })

  expect_true(all(bool_vec))
})

#########################################

## determine_minimium_neighborhood_size is correct

test_that("determine_minimium_neighborhood_size works", {
  mat <- cbind(1:10, 1:10)
  res <- determine_minimium_neighborhood_size(mat, verbose = F)

  expect_true(length(res) == 1)
  expect_true(res > 0)
  expect_true(res %% 1 == 0)
})

test_that("determine_minimium_neighborhood_size is correct", {
  set.seed(10)
  mat <- cbind(rep(0:1, each = 5), rep(0:1, each = 5))
  mat <- mat + rnorm(20, sd = 0.01)
  res <- determine_minimium_neighborhood_size(mat, verbose = F)

  expect_true(res == 5)
})

##############################

## .compute_pairwise_purity is correct

test_that(".compute_pairwise_purity works", {
  mat <- cbind(1:10, 1:10)
  neighborhood_size <- 3
  g <- .construct_neighborhood_graph(mat, neighborhood_size)
  cluster_labels <- rep(1:2, each = 5)

  res <- .compute_pairwise_purity(g, 1, 5, cluster_labels)

  expect_true(length(res) == 1)
  expect_true(res >= 0)
  expect_true(res <= 1)
})

test_that(".compute_pairwise_purity can get 0", {
  mat <- cbind(1:10, 1:10)
  neighborhood_size <- 3
  g <- .construct_neighborhood_graph(mat, neighborhood_size)
  cluster_labels <- c(1,1,1,1,2, 2,2,1,1,1)

  res <- .compute_pairwise_purity(g, 4, 8, cluster_labels)

  expect_true(res == 0)
})

######################################3

## compute_purity is correct

test_that("compute_purity works", {
  mat <- cbind(1:10, 1:10)
  neighborhood_size <- 3
  cluster_labels <- c(1,1,1,1,2, 2,2,1,1,1)

  res <- compute_purity(mat, cluster_labels, neighborhood_size)

  expect_true(is.list(res))
  expect_true(is.list(res$value_list))
  expect_true(length(res$avg_val) == 1)
})
