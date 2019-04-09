context("Test slingshot curves")

## .initialize_weight_matrix is correct

test_that(".initialize_weight_matrix works", {
  set.seed(10)
  cluster_labels <- sample(1:10, 200, replace = T)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
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
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
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
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
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
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
  pcurve_list <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list

  sample_idx <- .determine_idx_lineage(lineages[[1]], cluster_mat)
  res <- .smoother_func(pcurve_list[[1]]$lambda, dat[sample_idx,])

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == dim(dat[sample_idx,])))
})

########################

## .construct_average_curve is correct

test_that(".construct_average_curve works", {
  set.seed(10)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
  cluster_mat <- .construct_cluster_matrix(cluster_labels)
  k <- ncol(cluster_mat)
  centers <- .compute_cluster_center(dat, cluster_mat)
  W <- .initialize_weight_matrix(cluster_mat, lineages)
  cluster_vec <- 1:ncol(cluster_mat)
  s_list <- .initial_curve_fit(lineages, cluster_vec, centers)
  pcurve_list <- .refine_curve_fit(dat, s_list, lineages, W, cluster_mat)$pcurve_list

  res <- .construct_average_curve(pcurve_list, dat)

  expect_true(class(res) == "principal_curve")
  expect_true(all(!is.na(res$idx)))
})

###############################

## .refine_curve_fit is correct

test_that(".refine_curve_fit works", {
  set.seed(10)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 5), diag(5))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
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
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
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
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
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
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
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
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)
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

######################

## .get_curves is correct

test_that(".get_curves works", {
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
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)

  res <- .get_curves(dat, cluster_labels, lineages)

  expect_true(length(res) == length(lineages))
  expect_true(all(sapply(res, class) == "principal_curve"))

  #plot(dat[,1], dat[,2], asp = T)
  #lines(res$Curve1, col = "red", lwd = 2); lines(res$Curve2, col = "blue", lwd = 2)
})


test_that(".get_curves works for a harder example", {
  set.seed(20)
  cluster_labels <- rep(1:5, each = 20)
  dat <- MASS::mvrnorm(100, rep(0, 2), diag(2))
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)

  res <- .get_curves(dat, cluster_labels, lineages)

  expect_true(length(res) == length(lineages))
  expect_true(all(sapply(res, class) == "principal_curve"))

  #plot(dat[,1], dat[,2], asp = T)
  #for(i in 1:length(res)){lines(res[[i]], col = i, lwd = 2)}
})

test_that(".get_curves finds reasonable curves", {
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
  lineages <- .get_lineages(dat, cluster_labels, starting_cluster = 1)

  res <- .get_curves(dat, cluster_labels, lineages)

  bool_vec <- sapply(1:length(res), function(x){
    all(sapply(1:(length(lineages[[x]])-1), function(y){
      mean(res[[x]]$lambda[((y-1)*50+1):(y*50)]) < mean(res[[x]]$lambda[(y*50+1):((y+1)*50)])
    }))
  })

  expect_true(all(bool_vec))
})

##############################

## slingshot works

test_that("slingshot works", {
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
  res <- slingshot(dat, cluster_labels, starting_cluster = 1)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("lineages", "curves")))

  #plot(dat[,1], dat[,2], asp = T)
  #for(i in 1:length(res$curves)){lines(res$curves[[i]], col = i+1, lwd = 2)}
})


test_that("slingshot can give sensible lambdas", {
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
  res <- slingshot(dat, cluster_labels, starting_cluster = 1)

  bool_vec <- sapply(1:length(res$lineages), function(x){
    lin <- res$lineages[[x]]
    all(sapply(1:(length(lin)-1), function(y){
      idx1 <- which(cluster_labels == lin[y])
      idx2 <- which(cluster_labels == lin[y+1])

      mean(res$curves[[x]]$lambda[which(res$curves[[x]]$idx %in% idx1)]) <
        mean(res$curves[[x]]$lambda[which(res$curves[[x]]$idx %in% idx2)])
    }))
  })

  expect_true(all(bool_vec))
})


test_that("slingshot works with an artifical example", {
  set.seed(10)
  cluster_labels <- sample(1:10, 200, replace = T)
  dat <- MASS::mvrnorm(200, rep(0, 5), diag(5))

  res <- slingshot(dat, cluster_labels, starting_cluster = 1)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("lineages", "curves")))
})

############

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




