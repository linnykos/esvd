context("Test trajectory manipulation")

## .extract_pseudotimes is correct

test_that(".extract_pseudotimes works", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       60,80, 25,100,
                       40,10, 60,80,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_each <- 50
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
  }))
  cluster_labels <- rep(1:4, each = 50)
  slingshot_res <- slingshot(dat, cluster_labels, starting_cluster = 1)

  res <- .extract_pseudotimes(slingshot_res)

  expect_true(length(res) == 2)
  expect_true(is.list(res))
  expect_true(nrow(res[[1]]) == length(slingshot_res$curves[[1]]$idx))
  expect_true(nrow(res[[2]]) == length(slingshot_res$curves[[2]]$idx))
  expect_true(all(sapply(res, function(x){sort(colnames(x)) == sort(c("cell_idx", "pseudotime", "dist_to_curve"))})))
})

test_that(".extract_pseudotimes works with upsampling", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       60,80, 25,100,
                       40,10, 60,80,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_vec <- c(30,40,50,60)
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_vec[x])
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_vec[x], sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_vec[x], sd = 0.1))
  }))
  cluster_labels <- unlist(lapply(1:length(n_vec), function(i){rep(i, n_vec[i])}))
  slingshot_res <- slingshot(dat, cluster_labels, starting_cluster = 1, upscale_factor = 1)

  # check the cluster sizes
  tab <- table(cluster_labels[slingshot_res$idx])
  expect_true(tab[1]+tab[2] == tab[3])
  expect_true(tab[1]+tab[2] == tab[4])

  res <- .extract_pseudotimes(slingshot_res)

  expect_true(length(unique(res[[1]]$cell_idx)) == length(res[[1]]$cell_idx))
  expect_true(length(unique(res[[2]]$cell_idx)) == length(res[[2]]$cell_idx))
  expect_true(all(sort(res[[1]]$cell_idx) == c(1:(30+40+50))))
  expect_true(all(sort(res[[2]]$cell_idx) == c(1:(30+40), (30+40+50+1):(30+40+50+60))))
})

#########################

## .compile_common_cells is correct

test_that(".compile_common_cells works", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       60,80, 25,100,
                       40,10, 60,80,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_vec <- c(30,40,50,60)
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_vec[x])
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_vec[x], sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_vec[x], sd = 0.1))
  }))
  cluster_labels <- unlist(lapply(1:length(n_vec), function(i){rep(i, n_vec[i])}))
  slingshot_res <- slingshot(dat, cluster_labels, starting_cluster = 1, upscale_factor = 1)
  df_list <- .extract_pseudotimes(slingshot_res)

  res <- .compile_common_cells(df_list)

  expect_true(is.data.frame(res))
  expect_true(ncol(res) == 4)
  expect_true(all(sort(colnames(res)) == sort(c("cell_idx", "pseudotime", "dist_to_curve", "consensus"))))
  expect_true(all(res$cell_idx == 1:(30+40)))
})

##########################

## .compile_common_cells is correct

test_that(".compile_common_cells works", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       60,80, 25,100,
                       40,10, 60,80,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_vec <- c(30,40,50,60)
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_vec[x])
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_vec[x], sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_vec[x], sd = 0.1))
  }))
  cluster_labels <- unlist(lapply(1:length(n_vec), function(i){rep(i, n_vec[i])}))
  slingshot_res <- slingshot(dat, cluster_labels, starting_cluster = 1, upscale_factor = 1)
  df_list <- .extract_pseudotimes(slingshot_res)

  res <- .compile_common_cells(df_list)

  expect_true(is.data.frame(res))
  expect_true(ncol(res) == 4)
  expect_true(all(sort(colnames(res)) == sort(c("cell_idx", "pseudotime", "dist_to_curve", "consensus"))))
  expect_true(all(res$cell_idx == 1:(30+40)))
})

###############################

## .compile_unique_cells is correct

test_that(".compile_unique_cells works", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       60,80, 25,100,
                       40,10, 60,80,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_vec <- c(30,40,50,60)
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_vec[x])
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_vec[x], sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_vec[x], sd = 0.1))
  }))
  cluster_labels <- unlist(lapply(1:length(n_vec), function(i){rep(i, n_vec[i])}))
  slingshot_res <- slingshot(dat, cluster_labels, starting_cluster = 1, upscale_factor = 1)
  df_list <- .extract_pseudotimes(slingshot_res)

  res <- .compile_unique_cells(df_list)

  expect_true(is.list(res))
  expect_true(all(sapply(res, is.data.frame)))
  expect_true(all(sapply(res, ncol) == 3))
  expect_true(all(sapply(res, function(x){
    all(sort(colnames(x)) == sort(c("cell_idx", "pseudotime", "dist_to_curve")))
  })))
  expect_true(all(res[[1]]$cell_idx == c((30+40+1):(30+40+50))))
  expect_true(all(res[[2]]$cell_idx == c((30+40+50+1):(30+40+50+60))))
})

##################################3

## construct_pseudotime_trajectory_matrix is correct

test_that("construct_pseudotime_trajectory_matrix works", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       60,80, 25,100,
                       40,10, 60,80,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_vec <- c(30,40,50,60)
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_vec[x])
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_vec[x], sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_vec[x], sd = 0.1))
  }))
  cluster_labels <- unlist(lapply(1:length(n_vec), function(i){rep(i, n_vec[i])}))
  slingshot_res <- slingshot(dat, cluster_labels, starting_cluster = 1, upscale_factor = 1)

  res <- construct_pseudotime_trajectory_matrix(slingshot_res, cluster_labels)

  expect_true(is.data.frame(res))
  expect_true(all(res$cell_idx == c(1:nrow(dat))))
  expect_true(all(sort(colnames(res)) == sort(c("cell_idx", "pseudotime", "dist_to_curve", "consensus", "status", "cluster_labels"))))
})

#####################################

## partition_cells_using_pseudotime is correct

test_that("partition_cells_using_pseudotime works", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       60,80, 25,100,
                       40,10, 60,80,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_vec <- c(30,40,50,60)
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_vec[x])
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_vec[x], sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_vec[x], sd = 0.1))
  }))
  cluster_labels <- unlist(lapply(1:length(n_vec), function(i){rep(i, n_vec[i])}))
  slingshot_res <- slingshot(dat, cluster_labels, starting_cluster = 1, upscale_factor = 1)
  pseudotime_df <- construct_pseudotime_trajectory_matrix(slingshot_res, cluster_labels)

  res <- partition_cells_using_pseudotime(pseudotime_df, trajectory_1_clusters = 3,
                                          trajectory_2_clusters = 4)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == c("cell_idx_common", "cell_idx_traj1", "cell_idx_traj2")))
  expect_true(all(sort(res$cell_idx_common) %in% c(1:(30+40))))
  expect_true(all(sort(res$cell_idx_traj1) %in% c((30+40+1):(30+40+50))))
  expect_true(all(sort(res$cell_idx_traj2) %in% c((30+40+50+1):(30+40+50+60))))
})




