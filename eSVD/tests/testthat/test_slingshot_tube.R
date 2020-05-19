context("Test slingshot tube")

## .discretize_curve_by_pseudotime is correct

test_that(".discretize_curve_by_pseudotime works", {
  s_mat <- matrix(1:200, nrow = 40, ncol = 5)
  pseudotime_vec <- 1:40
  resolution <- 10

  res <- .discretize_curve_by_pseudotime(s_mat, pseudotime_vec, resolution)

  expect_true(length(res) == 2)
  expect_true(all(sort(names(res)) == c("lambda", "s")))
  expect_true(ncol(res$s) == ncol(s_mat))
  expect_true(all(res$s %in% s_mat))
  expect_true(length(res$lambda) == nrow(res$s))
  expect_true(all(res$lambda %in% pseudotime_vec))
})

test_that(".discretize_curve_by_pseudotime works if resolution is equal to nrow(s_mat)", {
  s_mat <- matrix(1:200, nrow = 40, ncol = 5)
  pseudotime_vec <- 1:40
  resolution <- 40

  res <- .discretize_curve_by_pseudotime(s_mat, pseudotime_vec, resolution)

  expect_true(length(res) == 2)
  expect_true(all(sort(names(res)) == c("lambda", "s")))
  expect_true(ncol(res$s) == ncol(s_mat))
  expect_true(all(res$s %in% s_mat))
  expect_true(length(res$lambda) == nrow(res$s))
  expect_true(all(res$lambda %in% pseudotime_vec))
})

test_that(".discretize_curve_by_pseudotime removes duplicates", {
  set.seed(10)
  s_mat <- matrix(1:200, nrow = 40, ncol = 5)
  pseudotime_vec <- cumsum(c(0, abs(rnorm(39))))
  resolution <- 10

  res <- .discretize_curve_by_pseudotime(s_mat, pseudotime_vec, resolution)

  s_mat2 <- s_mat
  for(i in 1:20){
    s_mat2 <- rbind(s_mat2[1,], s_mat2)
  }
  pseudotime_vec2 <- c(rep(0, 20), pseudotime_vec)
  res2 <- .discretize_curve_by_pseudotime(s_mat2, pseudotime_vec2, resolution)

  expect_true(sum(abs(res$s - res2$s)) <= 1e-6)
  expect_true(sum(abs(res$lambda - res2$lambda)) <= 1e-6)
})

##################################

## bootstrap_curves is correct

test_that("bootstrap_curves works", {
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

  lineages <- list(c(1,2,3), c(1,2,4))
  res <- bootstrap_curves(dat, cluster_labels, lineages = lineages, trials = 10)

  expect_true(length(res) == 10)
  expect_true(all(sapply(res, length) == 2))

  # for(i in 1:10){
  #   plot(dat[,1], dat[,2], asp = T, col = cluster_labels, pch = 16, main = i)
  #   for(j in 1:2){
  #     ord <- res[[i][[j]]$ord
  #     lines(res[[i]][[j]]$s[ord,1], res[[i]][[j]]$s[ord,2], lwd = 3, col = "white")
  #     lines(res[[i]][[j]]$s[ord,1], res[[i]][[j]]$s[ord,2])
  #   }
  # }
})

############################

## .compute_l2_curve is correct

test_that(".compute_l2_curve works", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       60,80, 25,100,
                       40,10, 60,80,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_each <- 25
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
  }))
  cluster_labels <- rep(1:4, each = n_each)

  set.seed(10)
  slingshot_res <- slingshot(dat, cluster_labels, starting_cluster = 1)
  lineages <- slingshot_res$lineages
  set.seed(10)
  trials <- 10
  bootstrap_res <- bootstrap_curves(dat, cluster_labels, lineages = lineages, trials = trials)

  mat <- slingshot_res$curves[[1]]$s
  mat_collection <- lapply(bootstrap_res, function(curves){curves[[1]]$s})

  res <- .compute_l2_curve(mat, mat_collection)

  expect_true(all(dim(res) == c(nrow(mat), trials)))
})

##############################

## compute_curve_sd is correct

test_that("compute_curve_sd works", {
  set.seed(10)
  cell_pop <- matrix(c(4,10, 25,100,
                       60,80, 25,100,
                       40,10, 60,80,
                       60,80, 100,25)/10, nrow = 4, ncol = 4, byrow = T)
  h <- nrow(cell_pop)
  n_each <- 25
  dat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = 0.1),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = 0.1))
  }))
  cluster_labels <- rep(1:4, each = n_each)

  set.seed(10)
  slingshot_res <- slingshot(dat, cluster_labels, starting_cluster = 1)
  lineages <- slingshot_res$lineages
  set.seed(10)
  trials <- 10
  bootstrap_res <- bootstrap_curves(dat, cluster_labels, lineages = lineages, trials = trials)

  res <- compute_curve_sd(slingshot_res$curves, bootstrap_res)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == c("mat_list", "target_curve_list")))
  expect_true(length(res$target_curve_list) == 2)
  expect_true(length(res$mat_list) == 2)
  expect_true(all(dim(res$mat_list[[1]]) == c(75, trials)))
  expect_true(all(dim(res$mat_list[[2]]) == c(75, trials)))
})

