context("Test slingshot lineages")

## .covariance_distance is correct

test_that(".covariance_distance works", {
  mean_vec1 <- c(0,0)
  cov_mat1 <- diag(c(1,2))
  n1 <- 10
  mean_vec2 <- c(1,0)
  cov_mat2 <- matrix(c(2,1,1,2), 2, 2)
  n2 <- 20

  res <- .covariance_distance(mean_vec1, cov_mat1, n1, mean_vec2, cov_mat2, n2)
  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
})

test_that(".covariance_distance is actually the test statistic for a hotelling test with unequal variances", {
  # results from http://www.real-statistics.com/multivariate-statistics/hotellings-t-square-statistic/hotellings-t-square-unequal-covariance-matrices/

  mean_vec1 <- c(0,0)
  mean_vec2 <- c(1,0)
  cov_mat1 <- diag(c(1,2))
  cov_mat2 <- matrix(c(2,1,1,2), 2, 2)
  n1 <- 100
  n2 <- 200
  trials <- 500

  test_vec1 <- sapply(1:trials, function(trial){
    set.seed(trial)
    dat1 <- MASS::mvrnorm(n1, mu = mean_vec1, Sigma = cov_mat1)
    dat2 <- MASS::mvrnorm(n2, mu = mean_vec1, Sigma = cov_mat2)
    .covariance_distance(colMeans(dat1), stats::cov(dat1), nrow(dat1),
                         colMeans(dat2), stats::cov(dat2), nrow(dat2))
  })

  test_vec2 <- sapply(1:trials, function(trial){
    set.seed(trial)
    dat1 <- MASS::mvrnorm(n1, mu = mean_vec1, Sigma = cov_mat1)
    dat2 <- MASS::mvrnorm(n2, mu = mean_vec2, Sigma = cov_mat2)
    .covariance_distance(colMeans(dat1), stats::cov(dat1), nrow(dat1),
                         colMeans(dat2), stats::cov(dat2), nrow(dat2))
  })

  set.seed(10)
  ideal_vec <- stats::rchisq(trials, df = length(mean_vec1))

  expect_true(sum(abs(sort(ideal_vec) - sort(test_vec1))) <= sum(abs(sort(ideal_vec) - sort(test_vec2))))
})

##########3#########

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

##################################

## .enumerate_dist_from_trees is correct

test_that(".enumerate_dist_from_trees works", {
  set.seed(10)
  dat <- MASS::mvrnorm(10, mu = c(0,0), Sigma = diag(2))
  dist_mat <- as.matrix(stats::dist(dat))
  tree_list <- list(c(1,2,5,10), c(1,2,5,8), c(1,2,4,6))

  res <- .enumerate_dist_from_trees(dist_mat, tree_list)

  expect_true(is.matrix(res))
  expect_true(ncol(res) == 3)
  expect_true(nrow(res) == sum(sapply(tree_list, function(tree){length(tree)-1})))

  for(i in 1:length(tree_list)){
    for(j in 1:(length(tree_list[[i]])-1)){
      node_pair <- tree_list[[i]][c(j,j+1)]
      idx <- intersect(which(res[,1] == node_pair[1]), which(res[,2] == node_pair[2]))
      expect_true(length(idx) >= 1)

      expect_true(abs(res[idx[1],3] - dist_mat[node_pair[1], node_pair[2]]) <= 1e-5)
    }
  }
})

###################

## .enumerate_dist_between_levels is correct

test_that(".enumerate_dist_between_levels works", {
  set.seed(15)
  dat <- MASS::mvrnorm(10, mu = c(0,0), Sigma = diag(2))
  dist_mat <- as.matrix(stats::dist(dat))
  tree_list <- list(c(1,2,5,10), c(1,2,5,8), c(1,2,4,6))
  cluster_vec <- c(3,7)

  res <- .enumerate_dist_between_levels(dist_mat, tree_list, cluster_vec)

  expect_true(is.matrix(res))
  expect_true(ncol(res) == 3)
  expect_true(nrow(res) == length(tree_list)*length(cluster_vec))

  for(i in 1:length(tree_list)){
    for(j in 1:length(cluster_vec)){
      node_pair <- c(tree_list[[i]][length(tree_list[[i]])], cluster_vec[j])
      idx <- intersect(which(res[,1] == node_pair[1]), which(res[,2] == node_pair[2]))
      expect_true(length(idx) == 1)

      expect_true(abs(res[idx[1],3] - dist_mat[node_pair[1], node_pair[2]]) <= 1e-5)
    }
  }
})

####################

## .enumerate_dist_within_levels is correct

test_that(".enumerate_dist_within_levels works", {
  set.seed(20)
  dat <- MASS::mvrnorm(10, mu = c(0,0), Sigma = diag(2))
  dist_mat <- as.matrix(stats::dist(dat))
  cluster_vec <- c(3,7,9)

  res <- .enumerate_dist_within_levels(dist_mat, cluster_vec)

  expect_true(is.matrix(res))
  expect_true(ncol(res) == 3)
  expect_true(nrow(res) == length(cluster_vec)*(length(cluster_vec)-1))

  for(i in 1:length(cluster_vec)){
    for(j in 1:length(cluster_vec)){
      if(cluster_vec[i] == cluster_vec[j]) next()

      node_pair <- cluster_vec[c(i,j)]
      idx <- intersect(which(res[,1] == node_pair[1]), which(res[,2] == node_pair[2]))
      expect_true(length(idx) == 1)

      expect_true(abs(res[idx[1],3] - dist_mat[node_pair[1], node_pair[2]]) <= 1e-5)
    }
  }
})

#############################

## .find_all_unique_paths is correct

test_that(".find_all_unique_paths works", {
  path_list <- list(c(1,3,5), c(1,3,5,8), c(1,2), c(1,2,10), c(1,2,10,14), c(1,3,7))

  res <- .find_all_unique_paths(path_list, starting_cluster = 1)

  expect_true(length(res) == 3)
  expect_true(any(sapply(res, function(x){length(x) == 4 && all(x == c(1,3,5,8))})))
  expect_true(any(sapply(res, function(x){length(x) == 4 && all(x == c(1,2,10,14))})))
  expect_true(any(sapply(res, function(x){length(x) == 3 && all(x == c(1,3,7))})))
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
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.5),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.5))
  }))
  cluster_labels <- rep(1:4, each = 50)

  ## plot(dat[,1], dat[,2], asp = T, pch = 16, col = cluster_labels)

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

