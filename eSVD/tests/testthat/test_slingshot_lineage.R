context("Test slingshot lineages")

## .compute_cluster_distances is correct

test_that(".compute_cluster_distances works", {
  set.seed(10)
  dat <- MASS::mvrnorm(n = 200, mu = rep(0, 5), Sigma = diag(5))
  cluster_labels <- rep(1:4, each = 50)

  res <- .compute_cluster_distances(dat, cluster_labels)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(4,4)))
})

##########

## .construct_graph is correct

test_that(".construct_graph works", {
  set.seed(10)
  dat <- MASS::mvrnorm(n = 200, mu = rep(0, 5), Sigma = diag(5))
  cluster_labels <- rep(1:4, each = 50)

  dist_mat <- .compute_cluster_distances(dat, cluster_labels)
  res <- .construct_graph(dist_mat)

  expect_true(class(res) == "igraph")
  expect_true(igraph::vcount(res) == 4)
  expect_true(igraph::ecount(res) == 4*3/2)
  expect_true(!igraph::is_directed(res))
})

test_that(".construct_graph has meaningful weights", {
  set.seed(20)
  dat <- MASS::mvrnorm(n = 200, mu = rep(0, 5), Sigma = diag(5))
  cluster_labels <- rep(1:4, each = 50)

  dist_mat <- .compute_cluster_distances(dat, cluster_labels)
  g <- .construct_graph(dist_mat)

  weight_mat <- cbind(igraph::as_edgelist(g, names = T), igraph::edge_attr(g, "weight", index = igraph::E(g)))

  combn_mat <- combn(4, 2)
  res <- lapply(1:ncol(combn_mat), function(x){
    i <- combn_mat[1,x]; j <- combn_mat[2,x]

    idx <- intersect(which(weight_mat[,1] == i), which(weight_mat[,2] == j))
    val1 <- abs(weight_mat[idx,3] - dist_mat[i,j])

    idx <- intersect(which(weight_mat[,2] == i), which(weight_mat[,1] == j))
    val2 <- abs(weight_mat[idx,3] - dist_mat[i,j])

    c(val1, val2)
  })

  expect_true(all(sum(unlist(res)) <= 1e-6))
})

##################

## .construct_lineage_from_hierarchy is correct

test_that(".construct_lineage_from_hierarchy works", {
  set.seed(10)
  cluster_labels <- sample(1:10, 200, replace = T)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))
  cluster_group_list <- list(1, c(2:5), c(6:8), c(9:10))

  res <- .get_lineages(dat, cluster_labels, cluster_group_list, starting_cluster = 1)

  expect_true(is.list(res))
})

test_that(".construct_lineage_from_hierarchy outputs a lineage that respects the cluster_group_list", {
  set.seed(10)
  cluster_labels <- sample(1:10, 200, replace = T)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))
  cluster_group_list <- list(1, c(2:5), c(6:8), c(9:10))
  cluster_encoding <- c(1, rep(2,4), rep(3,3), rep(4,2))

  res <- .get_lineages(dat, cluster_labels, cluster_group_list, starting_cluster = 1)

  # make sure all clusters are somewhere
  expect_true(length(unique(unlist(res))) == 10)
  for(i in 1:length(res)){
    tmp <- cluster_encoding[res[[i]]]
    expect_true(all(diff(tmp) >= 0))
  }
})

test_that(".construct_lineage_from_hierarchy has hierarchical that matters", {
  # construct a very simple example with spherical gaussians
  set.seed(10)
  cluster_labels <- rep(1:4, each = 100)
  dat <- rbind(MASS::mvrnorm(100, rep(0, 2), diag(2)),
               MASS::mvrnorm(100, rep(2, 2), diag(2)),
               MASS::mvrnorm(100, rep(-2, 2), diag(2)),
               MASS::mvrnorm(100, c(0.5, 0), diag(2)))
  cluster_group_list <- list(1, c(2:3), 4)

  res <- .get_lineages(dat, cluster_labels, cluster_group_list, starting_cluster = 1)
  res2 <- .get_lineages(dat, cluster_labels, starting_cluster = 1)

  # if cluster_group_list is provided, get the desired order
  if(length(res[[1]]) == 2){
    expect_true(all(res[[1]] == c(1,3)))
    expect_true(all(res[[2]] == c(1,2,4)))
  } else {
    expect_true(all(res[[2]] == c(1,3)))
    expect_true(all(res[[1]] == c(1,2,4)))
  }

  # if cluster_group_list is not provided, simply go for the shortest path lineage
  if(length(res2[[1]]) == 2){
    expect_true(all(res2[[1]] == c(1,3)))
    expect_true(all(res2[[2]] == c(1,4,2)))
  } else {
    expect_true(all(res2[[2]] == c(1,3)))
    expect_true(all(res2[[1]] == c(1,4,2)))
  }
})


##################

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

  ### construct the distance matrix
  dist_mat <- .compute_cluster_distances(dat, cluster_labels)

  ### construct the spt
  g <- .construct_graph(dist_mat)
  res <- .construct_lineages(g, starting_cluster = 1)

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
  res <- .get_lineages(dat, cluster_labels, starting_cluster = 1)

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
  res <- .get_lineages(dat, cluster_labels, starting_cluster = 1)

  expect_true(length(res) == 2)
  expect_true(all(res[[1]] == c(1,3,2)) | all(res[[2]] == c(1,3,2)))
  expect_true(all(res[[1]] == c(1,3,4)) | all(res[[2]] == c(1,3,4)))
  expect_true(any(res[[1]] != res[[2]]))
})

test_that(".get_lineages work for an artifical example", {
  set.seed(10)
  cluster_labels <- sample(1:10, 200, replace = T)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))

  res <- .get_lineages(dat, cluster_labels, starting_cluster = 1)

  expect_true(is.list(res))
})
