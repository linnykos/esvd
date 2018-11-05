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

####################

## .construct_spt is correct

test_that(".construct_spt works", {
  set.seed(10)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))
  knn_graph <- .construct_knn_graph(dat, 5)
  res <- .construct_spt(knn_graph, k = 5, starting_cluster = 1)

  expect_true(class(res) == "igraph")
})

test_that(".construct_spt finds the right graph for a specific configuration", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       40,10, 60,80,
                       60,80, 25,100,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_each <- 50
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
  }))
  cluster_labels <- rep(1:4, each = 50)

  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  centers <- .compute_cluster_center(dat, cluster_mat)
  dat_augment <- rbind(centers, dat)
  knn <- 1
  while(TRUE){
    knn_graph <- .construct_knn_graph(dat_augment, knn = knn)
    if(igraph::components(knn_graph)$no == 1) break()
    knn <- knn + 1
  }
  res <- .construct_spt(knn_graph, k = 4, starting_cluster = 1)
  res <- as.matrix(igraph::as_adjacency_matrix(res))

  expect_true(res[1,3] == 1)
  expect_true(res[3,2] == 1)
  expect_true(res[3,4] == 1)
  expect_true(sum(res) == 6)
})

##############

## .construct_lineages is correct

test_that(".construct_lineages finds the right lineage for a specific configuration", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       40,10, 60,80,
                       60,80, 25,100,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_each <- 50
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
  }))
  cluster_labels <- rep(1:4, each = 50)

  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  centers <- .compute_cluster_center(dat, cluster_mat)
  dat_augment <- rbind(centers, dat)
  knn <- 1
  while(TRUE){
    knn_graph <- .construct_knn_graph(dat_augment, knn = knn)
    if(igraph::components(knn_graph)$no == 1) break()
    knn <- knn + 1
  }
  spt_graph <- .construct_spt(knn_graph, k = 4, starting_cluster = 1)
  res <- .construct_lineages(spt_graph, starting_cluster = 1)

  expect_true(length(res) == 2)
  expect_true(all(res[[1]] == c(1,3,2)) | all(res[[2]] == c(1,3,2)))
  expect_true(all(res[[1]] == c(1,3,4)) | all(res[[2]] == c(1,3,4)))
  expect_true(any(res[[1]] != res[[2]]))
})

###################

## .get_lineages is correct

test_that(".get_lineages works", {
  set.seed(10)
  cluster_labels <- sample(1:10, 200, replace = T)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))
  res <- .get_lineages(dat, cluster_labels, knn = NA, starting_cluster = 1)

  expect_true(is.list(res))
})

test_that(".get_lineages finds the right lineage for a specific configuration", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       40,10, 60,80,
                       60,80, 25,100,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_each <- 50
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
  }))
  cluster_labels <- rep(1:4, each = 50)
  res <- .get_lineages(dat, cluster_labels, knn = NA, starting_cluster = 1)

  expect_true(length(res) == 2)
  expect_true(all(res[[1]] == c(1,3,2)) | all(res[[2]] == c(1,3,2)))
  expect_true(all(res[[1]] == c(1,3,4)) | all(res[[2]] == c(1,3,4)))
  expect_true(any(res[[1]] != res[[2]]))
})
